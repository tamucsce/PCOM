#pragma once

#include <boost/uuid/random_generator.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <string>
#include <string_view>
#include <type_traits>
#include <ostream>

namespace ligero::util::uuids {

struct random_generator {
    auto operator()() const { return boost::uuids::random_generator()(); }
};

struct uuid {
    std::string data;

    uuid() : data() { }
    uuid(random_generator rand) : data(boost::uuids::to_string(rand())) { }
    explicit uuid(const boost::uuids::uuid& id) : data(boost::uuids::to_string(id)) { }

    template <typename Str,
              typename = std::enable_if_t<std::is_same_v<std::decay_t<Str>, std::string>>>
    explicit uuid(Str&& str) : data(std::forward<Str>(str)) { }

    uuid(const uuid&)               = default;
    uuid(uuid&&) noexcept           = default;
    uuid& operator=(const uuid&)    = default;
    uuid& operator=(uuid&) noexcept = default;

    template <typename Str,
              typename = std::enable_if_t<std::is_same_v<std::decay_t<Str>, std::string>>>
    uuid& operator=(Str&& str) { data = std::forward<Str>(str); return *this; }

    bool operator==(const uuid& rhs) const noexcept { return data == rhs.data; }
    bool operator!=(const uuid& rhs) const noexcept { return !(*this == rhs); }
    bool operator< (const uuid& rhs) const noexcept { return data < rhs.data; }
    bool operator<=(const uuid& rhs) const noexcept { return data <= rhs.data; }

    const std::string& to_string() const { return data; }
    std::string_view to_string_view() const noexcept { return data; }
    
    decltype(auto) begin() noexcept { return data.begin(); }
    decltype(auto) cbegin() const noexcept { return data.cbegin(); }
    decltype(auto) end() noexcept { return data.end(); }
    decltype(auto) cend() const noexcept { return data.cend(); }

    auto size() const { return data.size(); }

    template <typename Archive>
    void serialize(Archive& ar, const unsigned int /* version */) { ar & data; }

    template <typename Archive, typename Iter>
    Iter __serialize(Archive ar, Iter buf) {
        return ar(buf, data);
    }
};

std::ostream& operator<<(std::ostream& os, const uuid& u) {
    os << u.data;
    return os;
}

}
