#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <tuple>

/*  Disable OpenMP support when compiling to WASM */
#if !defined(__EMSCRIPTEN__)
#include <omp.h>
#endif

namespace prelude {


// ==================== Pair ==================== //
constexpr auto fst = [](auto&& p) { return std::forward<decltype(p)>(p).first; };
constexpr auto snd = [](auto&& p) { return std::forward<decltype(p)>(p).second; };

// ==================== Tuple ==================== //
constexpr inline void tuple_head(std::tuple<>) { }
constexpr inline void tuple_tail(std::tuple<>) { }

template <typename Head, typename... Tail>
constexpr inline decltype(auto) tuple_head(const std::tuple<Head, Tail...>& tuple) {
    return std::get<Head>(tuple);
}

template <typename Head, typename... Tail>
constexpr inline decltype(auto) tuple_tail(std::tuple<Head, Tail...>&& tuple) {
    return std::apply([](auto&&, auto&&... tail) {
        return std::make_tuple<Tail...>(std::forward<Tail>(tail)...); }
        , std::forward<std::tuple<Head, Tail...>>(tuple));
}

// ==================== Category ==================== //
constexpr auto id = [](auto &&t) { return std::forward<decltype(t)>(t); };

constexpr auto compose = [](auto &&f, auto &&g) {
  return [f = std::forward<decltype(f)>(f),
          g = std::forward<decltype(g)>(g)](auto&&... a) {
             return f(g(std::forward<decltype(a)>(a)...));
         };
};

constexpr auto flip = [](auto&& f) {
    return [f = std::forward<decltype(f)>(f)](auto&& a, auto&& b) {
        return f(std::forward<decltype(b)>(b), std::forward<decltype(a)>(a));
    };
};

constexpr auto konst = [](auto&& k, auto&&...) { return std::forward<decltype(k)>(k); };

// ==================== Zips ==================== //
template <typename V, class BinFunc>
auto zip_with(BinFunc&& func, const V& a, const V& b) {
    auto min = std::min(a.size(), b.size());
    std::decay_t<V> tmp(min);
    if constexpr (std::is_fundamental_v<typename V::value_type>) {
        #pragma omp parallel for simd
        for (auto i = 0; i < min; i++) {
            tmp[i] = func(a[i], b[i]);
        }
    }
    else {
        #pragma omp parallel for
        for (auto i = 0; i < min; i++) {
            tmp[i] = func(a[i], b[i]);
        }
    }
    return tmp;
}

// operator+
// ------------------------------------------------------------
template <typename T = void>
struct plus : std::plus<T> {
    using std::plus<T>::operator();
};

template <typename T>
struct plus<std::pair<T, T>> : plus<T> {
    constexpr std::pair<T, T> operator()(const std::pair<T, T>& p0, const std::pair<T, T>& p1) {
        return { plus<T>::operator()(p0.first, p1.first), plus<T>::operator()(p0.second, p1.second) };
    }
};

template <typename T, typename V>
struct plus<std::pair<T, V>> : plus<T>, plus<V> {
    constexpr std::pair<T, V> operator()(const std::pair<T, V>& p0, const std::pair<T, V>& p1) {
        return { plus<T>::operator()(p0.first, p1.first), plus<V>::operator()(p0.second, p1.second) };
    }
};

template <template <typename...> typename C, typename... T>
struct plus<C<T...>> : std::plus<typename C<T...>::value_type> {
    constexpr C<T...> operator()(const C<T...>& c0, const C<T...>& c1) const {
        return prelude::zip_with(plus<typename C<T...>::value_type>(), c0, c1);
    }
};

template <>
struct plus<void> {
    template <typename T>
    constexpr T operator()(const T& x, const T& y) const {
        return prelude::plus<T>{}(x, y);
    }
};

template <typename T>
T operator+(const T& x, const T& y) {
    return prelude::plus<T>{}(x, y);
}

// operator^
// ------------------------------------------------------------
template <typename T = void>
struct bit_xor : std::bit_xor<T> {
    using value_type = T;
    using std::bit_xor<T>::operator();
};

template <typename T, typename V>
struct bit_xor<std::pair<T, V>> {
    using value_type = std::pair<T, V>;
    constexpr value_type operator()(const value_type& p0, const value_type& p1) const {
        return { bit_xor<T>()(p0.first, p1.first), bit_xor<V>()(p0.second, p1.second) };
    }
};

template <template <typename...> typename C, typename... T>
struct bit_xor<C<T...>> {
    using value_type = C<T...>;
    constexpr C<T...> operator()(const C<T...>& c0, const C<T...>& c1) const {
        return prelude::zip_with(bit_xor<typename C<T...>::value_type>(), c0, c1);
    }
};

template <>
struct bit_xor<void> {
    template <typename T>
    constexpr T operator()(const T& x, const T& y) const {
        return prelude::bit_xor<T>{}(x, y);
    }
};

constexpr auto zip_with_xor = [](const auto& x, const auto& y) {
    return zip_with(prelude::bit_xor<>(), x, y);
};

// Zip
// ------------------------------------------------------------
template <template <typename...> typename C,
          typename... T0,
          typename... T1,
          typename... T2>
auto zip3(C<T0...>&& c0, C<T1...>&& c1, C<T2...>&& c2) {
    assert(c0.size() == c1.size() && c1.size() == c2.size());
    auto iter0 = std::begin(c0), iter1 = std::begin(c1), iter2 = std::begin(c2);

    C<std::tuple<typename C<T0...>::value_type,
                 typename C<T1...>::value_type,
                 typename C<T2...>::value_type>> result;
    for (auto i = 0; i < c0.size(); i++) {
        result.emplace_back(std::move(*iter0), std::move(*iter1), std::move(*iter2));
        std::advance(iter0, 1);
        std::advance(iter1, 1);
        std::advance(iter2, 1);
    }
    return result;
}

// Reduction
// ------------------------------------------------------------
template <typename Container, typename Func>
typename Container::value_type reduce_with(Func&& f, Container&& in) {
    auto size = in.size();
    if (size == 1) { return in[0]; }
    
    auto middle = static_cast<decltype(size)>(size / 2);

    std::decay_t<Container> left(middle), right(middle);
    std::copy(std::make_move_iterator(in.begin()),
              std::make_move_iterator(in.begin() + middle),
              left.begin());
    if (size % 2) {
        std::copy(std::make_move_iterator(in.begin() + middle),
                  std::make_move_iterator(std::prev(in.end())),
                  right.begin());
        auto result = f(left, right);
        result.emplace_back(std::move(*std::prev(in.end())));
        return reduce_with(std::forward<Func>(f), std::move(result));
    }
    else {
        std::copy(std::make_move_iterator(in.begin() + middle),
                  std::make_move_iterator(in.end()),
                  right.begin());
        return reduce_with(std::forward<Func>(f), f(left, right));
    }
}

template <typename It, typename Func>
typename std::iterator_traits<It>::value_type
reduce_with(Func&& f, It begin, It end) {
    size_t len = std::distance(begin, end);
    assert(len > 0);
    if (len == 1) { return *begin; }

    size_t middle = len / 2;
    
    #pragma omp parallel for
    for (auto i = 0; i < middle; i++) {
        *std::next(begin, i) = f(*std::next(begin, i), *std::next(begin, i + middle));
    }

    if (len % 2) {
        *begin = f(*begin, *std::next(begin, len - 1));
    }
    return reduce_with(std::forward<Func>(f), begin, std::next(begin, middle));
}

// Concat
// ------------------------------------------------------------

// Move variant
template <template <typename...> typename Container,
          typename T, typename... CArgs, typename... TArgs>
Container<T, TArgs...> concat(Container<Container<T, TArgs...>, CArgs...>&& c) {
    if (c.empty()) { return Container<T, TArgs...>{}; }
    if (c.size() == 1) { return *std::begin(c); }

    // General case
    auto&& init = std::move(*std::begin(c));
    for (auto it = std::next(std::begin(c), 1); it != std::end(c); std::advance(it, 1)) {
        std::copy(std::make_move_iterator(std::begin(*it)),
                  std::make_move_iterator(std::end(*it)),
                  std::back_inserter(init));
    }
    return std::move(init);
}

// Copy variant
template <template <typename...> typename Container,
          typename T, typename... CArgs, typename... TArgs>
Container<T, TArgs...> concat(const Container<Container<T, TArgs...>, CArgs...>& c) {
    if (c.empty()) { return Container<T, TArgs...>{}; }
    if (c.size() == 1) { return *std::begin(c); }

    // General case
    auto init = *std::begin(c);
    for (auto it = std::next(std::begin(c), 1); it != std::end(c); std::advance(it, 1)) {
        std::copy(std::begin(*it), std::end(*it), std::back_inserter(init));
    }
    return init;
}

// Unconcat
// ------------------------------------------------------------
template <template <typename...> typename C, typename T, typename... TArgs>
auto unconcat(const C<T, TArgs...>& c, size_t n) {
    C<C<T, TArgs...>> out;
    for (auto it = std::begin(c); it != std::end(c); std::advance(it, n)) {
        out.emplace_back(it, std::next(it, n));
    }
    return out;
}

// Transpose
// ------------------------------------------------------------
template <typename T>
T transpose(const T& t) {
    if (t.size() == 0)
        return t;
    assert(t[0].size() != 0);

    T out;
    out.resize(t[0].size());
    for (auto& o : out) { o.resize(t.size()); }
    for (auto i = 0; i < t.size(); i++) {
        for (auto j = 0; j < t[0].size(); j++) {
            out[j][i] = t[i][j];
        }
    }
    return out;
}

// ==================== Abstract nonsense ==================== //
template <template <typename...> class K, typename... Ts>
struct apply_type {
    using type = K<Ts...>;
};

template <template <typename...> class K, typename V> struct apply_strip;
template <template <typename...> class K, template <typename...> class V, typename... Ts>
struct apply_strip<K, V<Ts...>> {
    using type = K<Ts...>;
};

// ------------------------------------------------------------

/// Helper function for std::visit
template <typename... Ts>
struct _overloaded : Ts... { using Ts::operator()...; };
template<class... Ts> _overloaded(Ts...) -> _overloaded<Ts...>;  // deduction guide

} // namespace prelude
