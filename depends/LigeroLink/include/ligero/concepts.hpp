#pragma once

#include <chrono>
#include <type_traits>
#include <functional>

namespace ligero::aya {

template <class T, class U>
concept same_as = std::is_same_v<T, U> && std::is_same_v<U, T>;

template <class From, class To>
concept convertible_to = std::is_convertible_v<From, To> &&
    requires(std::add_rvalue_reference_t<From> (&f)()) { static_cast<To>(f()); };

template <class Derived, class Base>
concept derived_from = std::is_base_of_v<Base, Derived> &&
    std::is_convertible_v<const volatile Derived*, const volatile Base*>;


template <class F, class... Args>
concept invocable = requires(F&& f, Args&&... args) {
    std::invoke(std::forward<F>(f), std::forward<Args>(args)...);
    /* not required to be equality preserving */
};

template <class Func> concept callable = requires(Func &&f) { f(); };


template <typename Transport, typename... Args>
concept has_pack = requires(Transport&& t, Args&&... args) {
    { t.pack(std::forward<Args>(args)...) } -> convertible_to<std::string>;
};

template <typename Transport, typename... Args>
concept has_unpack = requires(Transport&& t, const std::string& str, Args&&... args) {
    { t.template unpack<Args...>(str) } -> same_as<std::tuple<Args...>>;
};


template <typename Transport, typename... Args>
concept has_send = has_pack<Transport, Args...> &&
    requires(Transport&& t, Args&&... args) { { t.send(std::forward<Args>(args)...) }; };

template <typename Transport, typename... Args>
concept has_receive = has_unpack<Transport, Args...> &&
    requires(Transport&& t, std::chrono::nanoseconds timeout, Args&&... args) {
        { t.template receive<Args...>(timeout) } -> same_as<std::tuple<Args...>>;
    };

template <typename Transport, typename... Args>
concept duplex_transport = has_send<Transport, Args...> && has_receive<Transport, Args...>;

} // namespace ligero::aya

