#pragma once

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/utility.hpp>

#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>

#include <ligero/exception.hpp>
#include <ligero/util/boost/portable_binary_iarchive.hpp>
#include <ligero/util/boost/portable_binary_oarchive.hpp>
#include <ligero/util/timer.hpp>
#include <ligero/prelude.hpp>
#include <ligero/serialization/naive.hpp>

namespace ligero::net {

using namespace std::string_literals;

/**
 * Basic text archive using Boost serialization. Default flag
 * has been set to omit serialization header.
 *
 * Note: Typically this is the most portable serialize policy,
 *       but also the most heavy one (~2.7x in size compared
 *       with binary serialization).
 **/
struct boost_text {
    using iarchive_t = boost::archive::text_iarchive;
    using oarchive_t = boost::archive::text_oarchive;

    auto make_iarchive(std::stringstream& ss, unsigned int flags = boost::archive::no_header) const
        { return iarchive_t(ss, flags); }
    auto make_oarchive(std::stringstream& ss, unsigned int flags = boost::archive::no_header) const
        { return oarchive_t(ss, flags); }

protected:
    ~boost_text() = default;
};

/**
 * Basic text archive using Boost serialization. Default flag
 * has been set to omit serialization header.
 *
 * @note Boost binary archive is NOT portable. It serialize `std::string`
 *       into different size in native 64 bits build / 32 bits WASM build,
 *       which means you'll get a deserialization error.
 *
 *       Solution: Use `deflate_serializer` with `boost_text`
 *       or use `boost_portable_binary`
 **/
struct boost_binary {
    using iarchive_t = boost::archive::binary_iarchive;
    using oarchive_t = boost::archive::binary_oarchive;

    auto make_iarchive(std::stringstream& ss, unsigned int flags = boost::archive::no_header) const
        { return iarchive_t(ss, flags); }
    auto make_oarchive(std::stringstream& ss, unsigned int flags = boost::archive::no_header) const
        { return oarchive_t(ss, flags); }

protected:
    ~boost_binary() = default;
};

/**
 * Basic portable binary archive. Default flag has been set to
 * omit the serialization header and serialize as little endian.
 *
 * @note This is the best solution for now and should be set to
 *       default in every transport.
 **/
struct boost_portable_binary {
    using iarchive_t = boost::archive::portable_binary_iarchive;
    using oarchive_t = boost::archive::portable_binary_oarchive;

    auto make_iarchive(std::stringstream& ss, unsigned int flags = static_cast<unsigned int>(boost::archive::no_header) | static_cast<unsigned int>(boost::archive::endian_little)) const
        { return iarchive_t(ss, flags); }
    auto make_oarchive(std::stringstream& ss, unsigned int flags = static_cast<unsigned int>(boost::archive::no_header) | static_cast<unsigned int>(boost::archive::endian_little)) const
        { return oarchive_t(ss, flags); }
};


/**
 * Basic serializer templated with archive policy.
 *
 * @note The return value of unpack() is a tuple even if only one type was
 *       expected.
 *
 * @param ArchivePolicy  How to serialize data (text/binary/compressed text/...)
 *                       default: `boost_portable_binary`
 **/
template <typename ArchivePolicy = boost_portable_binary>
struct serializer : public ArchivePolicy {
    using oarchive_t = typename ArchivePolicy::oarchive_t;  //< generic output archive type
    using iarchive_t = typename ArchivePolicy::iarchive_t;  //< generic input archive type
    
    serializer() = default;
    serializer(size_t allocate_hint) { }  /*  Unused */
    
    /**
     * Pack a set of values into a single string
     *
     * @param val  Set of values
     * @return     Serialized string using given serialize policy
     **/
    template <typename... Args>
    auto pack(Args&&... val) {
        auto __t = util::make_timer("Aya", "Net", "Pack");
        std::stringstream ss;
        oarchive_t oa = this->make_oarchive(ss);

        // C++17 Fold expression
        // pushing values one by one from left to right
        ( (oa << val), ...);

        return ss.str();
    }

    /**
     * Unpack values from a string
     *
     * Usage: auto [a, b, c] = unpack<A, B, C>(msg);
     *
     * @param msg     A string-like message
     * @param Args... Unpack types
     **/
    template <typename... Args>
    std::tuple<Args...> unpack(std::span<const char> msg) const {
        auto __t = util::make_timer("Aya", "Net", "Unpack");
        std::stringstream ss;
        ss.write(msg.data(), msg.size());
        iarchive_t ia = this->make_iarchive(ss);
        
        try {
            return unpack_go<Args...>(ia);
        }
        catch (const std::exception& e) {
            LIGERO_THROW( serialize_error()
                          << throw_reason(e.what())
                          << throw_error_code(error_code::deserialize_failed) );
        }
    }

    template <typename... Args>
    void unpack(std::span<const char> msg, Args&... args) const {
        auto __t = util::make_timer("Aya", "Net", "Unpack");
        std::stringstream ss;
        ss.write(msg.data(), msg.size());
        iarchive_t ia = this->make_iarchive(ss);

        try {
            ( (ia >> args), ...);
        }
        catch (const std::exception& e) {
            LIGERO_THROW( serialize_error()
                          << throw_reason(e.what())
                          << throw_error_code(error_code::deserialize_failed) );
        }
    }

protected:
     /**
     * Recursively unpack value from an archive.
     * 
     * @param ia  An archive contains serialized string
     **/
    template <typename T, typename... Args>
    std::tuple<T, Args...> unpack_go(iarchive_t& ia) const {
        T val;
        ia >> val;
        if constexpr (sizeof...(Args) == 0) {
                return { std::move(val) };
            }
        else {
            return std::tuple_cat(std::make_tuple(std::move(val)), unpack_go<Args...>(ia));
        }
    }
};

// ------------------------------------------------------------

// Exciting!
struct naive_serializer {
    
    naive_serializer() = default;
    naive_serializer(size_t allocate_hint) { buffer_.reserve(allocate_hint); }
    
    template <typename... Args>
    auto pack(Args&&... args) {
        auto __t = util::make_timer("Aya", "Net", "Pack");
        buffer_.clear();
        serialize(std::back_inserter(buffer_), std::forward<Args>(args)...);
        return std::span<const char>{ buffer_.data(), buffer_.size() };
    }

    template <typename... Args>
    auto unpack(std::span<const char> buf, Args&&... args) const {
        auto __t = util::make_timer("Aya", "Net", "Unpack");
        return deserialize(buf.data(), std::forward<Args>(args)...);
    }

    template <typename T, typename... Args>
    std::tuple<T, Args...> unpack(std::span<const char> buf) const {
        auto __t = util::make_timer("Aya", "Net", "Unpack");
        T val;
        auto it = deserialize(buf.data(), val);
        if constexpr (sizeof...(Args) == 0) {
                return { std::move(val) };
            }
        else {
            return std::tuple_cat(std::make_tuple(std::move(val)),
                                  unpack<Args...>(buf.subspan(std::distance(buf.data(), it))));
        }
    }

protected:
    std::string buffer_;
};

} // namespace ligero::aya::net

