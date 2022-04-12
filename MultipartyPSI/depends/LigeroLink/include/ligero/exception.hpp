#pragma once

#include <boost/exception/all.hpp>
#include <cstdint>
#include <exception>

#define LIGERO_THROW(e) BOOST_THROW_EXCEPTION(e)

namespace ligero {

// --------------------  Exception Base  -------------------- //
struct exception       : virtual std::exception, virtual boost::exception { };
struct runtime_error   : virtual exception { };

// --------------------     IO Error     -------------------- //
struct system_error    : virtual runtime_error { };
struct io_error        : virtual runtime_error { };
struct input_error     : virtual io_error { };
struct network_error   : virtual io_error { };
struct serialize_error : virtual io_error { };

// -------------------- High-level Error -------------------- //
struct logic_error     : virtual runtime_error { };
struct protocol_error  : virtual runtime_error { };
struct timeout_error   : virtual protocol_error, virtual io_error { };

// --------------------   Misc Error    -------------------- //
struct math_error      : virtual runtime_error { };
struct unknown_error   : virtual runtime_error { };


// -------------------- Wrapper Functions -------------------- //
template <class ErrorInfo, class E>
decltype(auto) get_info(E& e) {
    return ::boost::get_error_info<ErrorInfo, E>(e);
}

template <typename Exception>
decltype(auto) diagnostic_info(const Exception& e, bool verbose=true) {
    return ::boost::diagnostic_information(e, verbose);
}

// --------------------    Error Code    -------------------- //
enum class error_code  : uint16_t {
    unknown_error = 0,
    websocket_not_supported,
    socket_creation_failed,
    socket_open_failed,
    peer_kicked_out,
    peer_uuid_collision,
    invalid_input,
    file_not_found,
    unrecognized_file_format,
    receive_timeout,
    send_timeout,
    no_message_received,
    deserialize_failed,
    compress_failed,
    decompress_failed,
    use_before_initialize,
    out_of_sync,
    too_few_clients,
    duplicate_uuid,
    zmq_receive_error,
    zmq_send_error,
};

// -------------------- Additional Info -------------------- //
using throw_api_function = boost::error_info<struct throw_api_function_, std::string>;
using throw_error_code = boost::error_info<struct throw_error_code_, error_code>;
using throw_reason = boost::error_info<struct throw_reason_, std::string>;


void append_api_function(exception& e, std::string func) {
    if (auto* ptr = get_info<throw_api_function>(e)) {
        ptr->append(std::string(" <- ") + func);
    }
    else {
        e << throw_api_function(std::move(func));
    }
}

} // namespace ligero

