#pragma once

#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

#define TRACE BOOST_LOG_TRIVIAL(trace)
#define DEBUG BOOST_LOG_TRIVIAL(debug)
#define INFO BOOST_LOG_TRIVIAL(info)
#define WARNING BOOST_LOG_TRIVIAL(warning)
#define ERROR BOOST_LOG_TRIVIAL(error)
#define FATAL BOOST_LOG_TRIVIAL(fatal)

namespace logging = boost::log;

enum class log_level : unsigned char {
    disabled,
    debug_only,
    info_only,
    full,
};

void enable_logging() {
    if (!logging::core::get()->get_logging_enabled())
        logging::core::get()->set_logging_enabled(true);
}

void disable_logging() {
    if (logging::core::get()->get_logging_enabled())
        logging::core::get()->set_logging_enabled(false);
}

void set_logging_level(log_level level) {
    switch (level) {
    case log_level::disabled:
        disable_logging();
        break;
    case log_level::debug_only:
        enable_logging();
        logging::core::get()->set_filter
            (logging::trivial::severity == logging::trivial::debug);
        break;
    case log_level::info_only:
        enable_logging();
        logging::core::get()->set_filter
            (logging::trivial::severity == logging::trivial::info);
        break;
    case log_level::full:
        enable_logging();
        boost::log::core::get()->set_filter
            (logging::trivial::severity >= logging::trivial::trace);
        break;
    }
}
