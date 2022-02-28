#pragma once

#include <chrono>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <string_view>
#include <thread>
#include <vector>
#include <iostream>

#include <ligero/util/color.hpp>
#include <ligero/util/log.hpp>
#include <ligero/util/tree.hpp>
#include <ligero/util/singleton.hpp>

namespace ligero::util {

struct nano {
    using resolution = std::chrono::nanoseconds;
    using duration_t = std::chrono::duration<size_t, std::nano>;
    constexpr static auto time_unit = "ns";
};

struct milli {
    using resolution = std::chrono::milliseconds;
    using duration_t = std::chrono::duration<size_t, std::milli>;
    constexpr static auto time_unit = "ms";
};

template <class Resolution = milli>
struct high_resolution_timer {
    using timepoint_t = std::chrono::high_resolution_clock::time_point;
    using duration_t = typename Resolution::duration_t;
    
    auto now() const { return std::chrono::high_resolution_clock::now(); }
    auto duration(const timepoint_t& start, const timepoint_t& end) const {
        return std::chrono::duration_cast<typename Resolution::resolution>(end - start);
    }
};


template <typename T, typename Name, typename... Args>
auto find_timer(tree<T>& current, Name&& name, Args&&... args) {
    auto it = std::find_if(current.children().begin(), current.children().end(),
                           [&name](auto&& child) {
                               return child.value().first == name;
                           });
        
    // not found - create new node
    if (it == current.children().end()) {
        it = current.emplace({name, {}});
    }
        
    if constexpr (sizeof...(args) == 0) {
        return it;
    }
    else {
        return find_timer(*it, std::forward<Args>(args)...);
    }
}

template <class Resolution>
struct basic_timer
{
    using lock_t = std::mutex;
    using lock_guard_t = std::scoped_lock<std::mutex>;
    using timepoint_t = typename high_resolution_timer<Resolution>::timepoint_t;
    using duration_seq = std::vector<size_t>;
    using timer_rep = std::pair<std::string_view, duration_seq>;

    using sumtree_node_t = std::tuple<std::string_view, size_t, size_t, size_t, size_t, size_t>;  // <name, sum_of_current, sum_of_child, min, max, avg>

    template <typename Name, typename... Args>
    struct timer_guard : public high_resolution_timer<Resolution> {
        timer_guard(tree<timer_rep>& root, lock_t& lock, Name&& name, Args&&... args)
            : stopped_(false),
              start_time_(this->now()),
              root_(root),
              args_(std::make_tuple(std::move(name), std::move(args)...))
            , lock_(lock) { }
        
        ~timer_guard() { stop(); }
        
        void stop() {
            if (stopped_) return;
            
            [[maybe_unused]] lock_guard_t lock_guard(lock_);
            auto lap = this->duration(start_time_, this->now());

            // Note: we have get latest reference every time we want
            //       to insert, because allocation might invalidate
            //       existing references.
            auto it = std::apply(
                [this](auto&&... args){
                    return find_timer(root_, std::forward<decltype(args)>(args)...);
                }, args_);
            it->value().second.emplace_back(lap.count());
            
            stopped_ = true;
        }
        
    private:
        bool stopped_;
        timepoint_t start_time_;
        tree<timer_rep>& root_;
        std::tuple<Name, Args...> args_;
        lock_t& lock_;
    };

    template <typename... Args>
    auto make_timer(Args... args) {
        return timer_guard<Args...>(timers_, mutex_, std::forward<Args>(args)...);
    }

    auto sum_of_timer(const tree<timer_rep>& root) const {
        auto go =
            [](const tree<timer_rep>& node,
               tree<sumtree_node_t>& sum,
               auto&& self)
                {
                    bool seq_empty = false;
                    size_t sum_of_current = 0;
                    size_t min = 0, max = 0, avg = 0;

                    if (node) {
                        auto& seq = node.value().second;
                        seq_empty = seq.empty();
                        sum_of_current = std::reduce(seq.begin(), seq.end(), sum_of_current);
                        if (!seq.empty()) {
                            min = *std::min_element(seq.begin(), seq.end());
                            max = *std::max_element(seq.begin(), seq.end());
                            avg = sum_of_current / seq.size();
                        }
                    }
                    // base case: node doesn't have child
                    if (!node.has_child()) {
                        sum.data() = node
                            ? std::make_optional(
                                std::make_tuple(node.value().first,
                                                sum_of_current,
                                                0,
                                                min, max, avg))
                            : std::nullopt;
                        return sum_of_current;
                    }
                    // has child node: run DFS
                    else {
                        size_t sum_of_children = 0;
                        for (auto& child : node.children()) {
                            auto sum_child = sum.emplace({"", 0, 0, 0, 0, 0});
                            sum_of_children += self(child, *sum_child, self);
                        }

                        if (node) {
                            sum.data() = {
                                node.value().first, sum_of_current, sum_of_children,
                                min, max, avg
                            };
                        }
                        return seq_empty ? sum_of_children : sum_of_current;
                    }
                };
        tree<sumtree_node_t> sumtree;
        go(root, sumtree, go);
        return sumtree;
    }

    template <typename Stream>
    Stream& __format(Stream& os,
                     const tree<sumtree_node_t>& st,
                     std::vector<std::string>& prefix) const
    {
        if (st) {
            for (size_t i = 0; i < prefix.size(); i++) { os << "    "; }
            for (auto pref : prefix) { os << pref << "."; }
        
            os << std::get<0>(st.value()) << ": ";

            auto [name, self, child, min, max, avg] = st.value();

            // Show current and children timer:
            // 1. current  = 0: show sum of child timer
            // 2. current != 0: show current timer
            if (self) {
                os << ANSI_GREEN << self << Resolution::time_unit << ANSI_RESET
                   << "    (min: "
                   << ANSI_GREEN << min << ANSI_RESET
                   << ", max: "
                   << ANSI_GREEN << max << ANSI_RESET
                   << ", avg: "
                   << ANSI_GREEN << avg << ANSI_RESET
                   << ")";
            }
            else {
                os << ANSI_GREEN << child << Resolution::time_unit << ANSI_RESET;
            }
            os << std::endl;

            prefix.emplace_back(name);
        }

        for (const auto& child : st.children()) { __format(os, child, prefix); }
        if (st && !prefix.empty()) { prefix.pop_back(); }
        return os;
    }

    template <typename Stream>
    Stream& format(Stream& os, const tree<timer_rep>& t) const {
        std::vector<std::string> prefix;
        auto sumtree = sum_of_timer(t);
        __format(os, sumtree, prefix);
        return os;
    }

    void print() const {
        std::cout << std::endl
                  << "============= Timing Info =============" << std::endl;
        
        format(std::cout, timers_);
        std::cout << "=======================================" << std::endl;
    }

public:
    tree<timer_rep> timers_;
    std::mutex mutex_;
};

using timer_t = singleton<basic_timer<milli>>;

template <typename... Args>
auto make_timer(Args... args) {
    return timer_t::instance().make_timer(std::forward<Args>(args)...);
}

void show_timer() { return timer_t::instance().print(); }

}  // namespace ligero::util
