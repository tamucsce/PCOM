#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdlib>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include <ligero/concepts.hpp>
#include <ligero/util/log.hpp>

namespace ligero::aya {

/*  Comply with Boost serialization's interface, actually unused  */
constexpr unsigned int naive_serialize_version = 100;

// Forward declaration
// ------------------------------------------------------------
template <typename Iter, typename T> struct do_serialize {
    static_assert(T::_, "Object is not serializable");
};

template <typename Iter, typename T> struct do_deserialize {
    static_assert(T::_, "Object is not serializable");
};

template <typename> struct proxy_serialize;
template <typename> struct proxy_deserialize;


// Concepts
// ------------------------------------------------------------
template <typename T, typename Iter>
concept Serializable = requires(T t, proxy_serialize<Iter> p, proxy_deserialize<Iter> dp, unsigned int v) {
    { t.serialize(p, v) } -> same_as<void>;
    { t.serialize(dp, v) } -> same_as<void>;
};

template <typename T>
concept IsFundamental = std::is_fundamental_v<T>;

template <typename T>
concept IsEnum = std::is_enum_v<T>;

template <typename T>
concept IsArray = std::is_array_v<T>;


// Proxy functors
// ------------------------------------------------------------
template <typename Iter>
struct proxy_serialize {
    proxy_serialize(Iter it) : it_(it) { }

    auto iter() { return it_; }
    
    template <typename... Args>
    constexpr Iter operator()(Args&&... args) {
        ( (it_ = do_serialize<Iter, std::remove_cv_t<std::remove_reference_t<Args>>>{}(
               it_, std::forward<Args>(args))), ... );
        return it_;
    }
    
protected:
    Iter it_;
};

template <typename Iter>
struct proxy_deserialize {
    proxy_deserialize(Iter it) : it_(it) { }

    auto iter() { return it_; }
    
    template <typename... Args>
    constexpr Iter operator()(Args&&... args) {
        ( (it_ = do_deserialize<Iter, std::remove_cv_t<std::remove_reference_t<Args>>>{}(
               it_, std::forward<Args>(args))), ... );
        return it_;
    }

protected:
    Iter it_;
};


template <typename Iter, typename T>
proxy_serialize<Iter>& operator&(proxy_serialize<Iter>& ar, const T& data) {
    ar(data);
    return ar;
}

template <typename Iter, typename T>
proxy_deserialize<Iter>& operator&(proxy_deserialize<Iter>& ar, T& data) {
    ar(data);
    return ar;
}


// Generic interfaces
// ------------------------------------------------------------
template <typename Iter, typename... Args>
Iter serialize(Iter buf, Args&&... args) {
    return proxy_serialize{buf}(std::forward<Args>(args)...);
}

template <typename Iter, typename... Args>
Iter deserialize(Iter buf, Args&&... args) {
    return proxy_deserialize{buf}(std::forward<Args>(args)...);
}


// Serialization
// ------------------------------------------------------------

// Primitive
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsFundamental<T>
struct do_serialize<Iter, T> {
    constexpr Iter operator()(Iter buffer, const T& data) {
        auto ptr = reinterpret_cast<const unsigned char *>(&data);
        return std::copy(ptr, ptr + sizeof(T), buffer);
    }
};

// Enum
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsEnum<T>
struct do_serialize<Iter, T> {
    constexpr Iter operator()(Iter buffer, const T& data) {
        return do_serialize<Iter, std::underlying_type_t<T>>{}(buffer, static_cast<std::underlying_type_t<T>>(data));
    }
};

// String
// ------------------------------------------------------------
template <typename Iter>
struct do_serialize<Iter, std::string> {
    constexpr Iter operator()(Iter buffer, const std::string& data) {
        auto it = do_serialize<Iter, uint64_t>{}(buffer, static_cast<uint64_t>(data.size()));
        return std::copy(data.begin(), data.end(), it);
    }   
};

// Custom
// ------------------------------------------------------------
template <typename Iter, typename T> requires Serializable<T, Iter>
struct do_serialize<Iter, T> {
    constexpr Iter operator()(Iter buffer, const T& data) {
        auto archive = proxy_serialize{buffer};
        const_cast<T&>(data).serialize(archive, naive_serialize_version);
        return archive.iter();
    }
};

// Vector
// ------------------------------------------------------------
template <typename Iter, typename T>
struct do_serialize<Iter, std::vector<T>> {
    constexpr Iter operator()(Iter buffer, const std::vector<T>& data) {
        uint64_t vec_len = data.size();
        auto it = do_serialize<Iter, uint64_t>{}(buffer, vec_len);

        for (auto i = 0; i < vec_len; i++) {
            it = do_serialize<Iter, T>{}(it, data[i]);
        }
        return it;
    }
};

// Pair
// ------------------------------------------------------------
template <typename Iter, typename T, typename V>
struct do_serialize<Iter, std::pair<T, V>> {
    constexpr Iter operator()(Iter buffer, const std::pair<T, V>& data) {
        auto it = do_serialize<Iter, T>{}(buffer, data.first);
        return do_serialize<Iter, V>{}(it, data.second);
    }
};

// Array
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsArray<T>
struct do_serialize<Iter, T> {
    constexpr Iter operator()(Iter buf, const T& data) {
        for (auto i = std::begin(data); i != std::end(data); i++) {
            buf = do_serialize<Iter, std::remove_cv_t<std::remove_reference_t<decltype(*i)>>>{}(buf, *i);
        }
        return buf;
    }
};



// Deserialization
// ------------------------------------------------------------

// Primitive
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsFundamental<T>
struct do_deserialize<Iter, T> {
    constexpr Iter operator()(Iter buf, T& out) {
        auto ptr = reinterpret_cast<unsigned char*>(&out);
        std::copy(buf, buf + sizeof(T), ptr);
        return buf + sizeof(T);
    }
};

// Enum
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsEnum<T>
struct do_deserialize<Iter, T> {
    constexpr Iter operator()(Iter buf, T& out) {
        return do_deserialize<Iter, std::underlying_type_t<T>>{}(buf, reinterpret_cast<std::underlying_type_t<T>&>(out));
    }
};

// String
// ------------------------------------------------------------
template <typename Iter>
struct do_deserialize<Iter, std::string> {
    constexpr Iter operator()(Iter buf, std::string& out) {
        uint64_t len = 0;
        auto it = do_deserialize<Iter, uint64_t>{}(buf, len);
        out.resize(len);
        std::copy(it, it + len, out.begin());
        return it + len;
    }
};

// Custom
// ------------------------------------------------------------
template <typename Iter, typename T> requires Serializable<T, Iter>
struct do_deserialize<Iter, T> {
    constexpr Iter operator()(Iter buf, T& out) {
        auto archive = proxy_deserialize{buf};
        out.serialize(archive, naive_serialize_version);
        return archive.iter();
    }
};

// Vector
// ------------------------------------------------------------
template <typename Iter, typename T>
struct do_deserialize<Iter, std::vector<T>> {
    constexpr Iter operator()(Iter buf, std::vector<T>& out) {
        uint64_t len = 0;
        auto it = do_deserialize<Iter, uint64_t>{}(buf, len);
        out.resize(len);

        for (auto i = 0; i < len; i++) {
            it = do_deserialize<Iter, T>{}(it, out[i]);
        }
        return it;
    }
};

// Pair
// ------------------------------------------------------------
template <typename Iter, typename T, typename V>
struct do_deserialize<Iter, std::pair<T, V>> {
    constexpr Iter operator()(Iter buf, std::pair<T, V>& out) {
        auto it = do_deserialize<Iter, T>{}(buf, out.first);
        return do_deserialize<Iter, V>{}(it, out.second);
    }
};

// Array
// ------------------------------------------------------------
template <typename Iter, typename T> requires IsArray<T>
struct do_deserialize<Iter, T> {
    constexpr Iter operator()(Iter buf, T& out) {
        for (auto i = std::begin(out); i != std::end(out); i++) {
            buf = do_deserialize<Iter, std::remove_cv_t<std::remove_reference_t<decltype(*i)>>>{}(buf, *i);
        }
        return buf;
    }
};

} // namespace ligero::aya

