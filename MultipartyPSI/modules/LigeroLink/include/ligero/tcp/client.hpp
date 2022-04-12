#pragma once

#include <chrono>
#include <cstdlib>
#include <string>

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>
#include <boost/system/error_code.hpp>
#include <boost/system/system_error.hpp>

#include <ligero/concepts.hpp>
#include <ligero/exception.hpp>
#include <ligero/prelude.hpp>
#include <ligero/util/log.hpp>
#include <ligero/util/timer.hpp>
#include <ligero/util/uuid.hpp>

namespace ligero::net {

using namespace boost;
using namespace std::chrono;
using namespace std::chrono_literals;

using asio::awaitable;
using asio::co_spawn;
using asio::detached;
using asio::use_awaitable;
using asio::ip::tcp;

/************************************************************
 * Basic class for general websocket client
 ************************************************************/
template <typename SerializePolicy>
struct basic_tcp_client {
    
    using socket_t = beast::tcp_stream;
    using message_t = std::string;
    using identity_t = util::uuids::uuid;

protected:
    size_t num_threads_ = 1;
    SerializePolicy archive_;       /**< Serialize archive  */
    identity_t identity_;           /**< Client's identity   */
    asio::io_context& executor_;    /**< IO context          */
    socket_t socket_;            /**< Socket              */
    beast::flat_buffer buffer_;     /**< Message buffer   */

public:
    basic_tcp_client(asio::io_context& exe)
        : basic_tcp_client(exe, "", 0) { }
    basic_tcp_client(asio::io_context& exe, std::string id)
        : basic_tcp_client(exe, std::move(id), 0) { }
    basic_tcp_client(asio::io_context& exe, std::string id, size_t allocate_hint)
        : archive_(allocate_hint)
        , identity_(id.empty() ? identity_t(util::uuids::random_generator()) : identity_t(id))
        , executor_(exe)
        , socket_(exe)
        { }

    virtual ~basic_tcp_client() = default;

    auto identity() const noexcept { return identity_; }
    void identity(identity_t id) { identity_ = id; }

    size_t num_threads() const noexcept { return num_threads_; }
    void num_threads(size_t num) { num_threads_ = num; }

    void launch_coroutine() {
        if (num_threads_ == 0) {
            ERROR << "Failed to get hardware concurrency, set it to 1";
            num_threads_ = 1;
        }

        asio::thread_pool pool(num_threads_);
        for (size_t i = 0; i < num_threads_; i++) {
            asio::post(pool, [this] { executor_.run(); });
        }
        
        pool.join();
        executor_.restart();
    }

    /**
     * Connect to remote server. Required before other operations
     *
     * @param connect_timeout  Connection timeout
     **/
    void connect(std::string endpoint, std::string port) {
        std::exception_ptr err;
        co_spawn(executor_,
                 coro_connect(endpoint, port),
                 [&err](std::exception_ptr eptr, int rc) {
                     err = eptr;
                 });

        launch_coroutine();

        if (err) {
            std::rethrow_exception(err);
        }
    }

    /**
     * Send message to remote
     *
     * @param ...  Data
     **/
    template <typename... Args>
    void send(Args&&... args) {
        auto view = archive_.pack(std::forward<Args>(args)...);
        auto buffer = asio::buffer(view.data(), view.size());

        uint64_t header = view.size();
        auto header_buffer = asio::buffer(
            reinterpret_cast<const char*>(&header), sizeof(uint64_t));
        
        std::exception_ptr err;
        co_spawn(executor_,
                 [&]() -> awaitable<void> {
                     co_await asio::async_write(socket_, header_buffer, use_awaitable);
                     co_await asio::async_write(socket_, buffer, use_awaitable);
                 },
                 [&err](std::exception_ptr e) { err = e; });

        launch_coroutine();

        if (err) {
            std::rethrow_exception(err);
        }
    }

    /**
     * Receive a single message from remote
     *
     * @param receive_timeout  Receive timeout
     * @return  Deserialized message
     **/
    template <typename... Args>
    std::tuple<Args...> receive(nanoseconds receive_timeout) {
        raw_receive(receive_timeout);
        auto t = archive_.template unpack<Args...>(
            { static_cast<const char *>(buffer_.data().data()), buffer_.data().size() });
        buffer_.consume(buffer_.size());
        return t;
    }
    

protected:
    /**
     * Raw receive without deserialize the message
     *
     * @param receive_timeout  Receive timeout
     * @return  A string like object
     **/
    void raw_receive(nanoseconds receive_timeout) {
        std::exception_ptr err;
        co_spawn(executor_,
                 [&]() -> awaitable<void> {
                     constexpr auto header_size = 8;  // uint64_t

                     buffer_.consume(buffer_.size());
                     
                     // Set timeout
                     beast::get_lowest_layer(socket_).expires_after(receive_timeout);

                     auto header_buf = buffer_.prepare(header_size);
                     auto n = co_await asio::async_read(socket_, header_buf, use_awaitable);
                     buffer_.commit(n);
                     
                     uint64_t message_size = 0;
                     std::copy(reinterpret_cast<const char*>(buffer_.data().data()),
                               reinterpret_cast<const char*>(buffer_.data().data()) + buffer_.data().size(),
                               reinterpret_cast<char*>(&message_size));
                     buffer_.consume(buffer_.size());

                     auto buf = buffer_.prepare(message_size);
                     n = co_await asio::async_read(socket_, buf, use_awaitable);
                     buffer_.commit(n);

                     // Reset timeout
                     beast::get_lowest_layer(socket_).expires_never();
                 },
                 [&err](std::exception_ptr e) { err = e; });
        launch_coroutine();

        if (err) {
            std::rethrow_exception(err);
        }
    }
    
    awaitable<int> coro_connect(std::string host, std::string port) {
        auto executor = co_await asio::this_coro::executor;

        tcp::resolver resolver(executor);
        const auto remote = co_await resolver.async_resolve(host, port, use_awaitable);

        std::string remote_url;
        for (;;) {
            try {
                auto ep = co_await asio::async_connect(beast::get_lowest_layer(socket_).socket(), remote, use_awaitable);
                remote_url = host + ":" + std::to_string(ep.port());
                DEBUG << "Connected to " << remote_url;
                break;
            }
            catch (system::system_error& e) {
                if (e.code() != system::errc::connection_refused) {
                    DEBUG << e.what();
                    LIGERO_THROW( system_error()
                                  << throw_reason(e.what())
                                  << throw_api_function("coro_connect") );
                }
            }


            asio::high_resolution_timer timer(socket_.get_executor());
            timer.expires_after(20ms);
            co_await timer.async_wait(use_awaitable);
        }

        {
            auto& tcp_socket = beast::get_lowest_layer(socket_).socket();
            tcp_socket.set_option(asio::socket_base::keep_alive{true});
            tcp_socket.set_option(tcp::no_delay{true});
        }

        co_return 0;
    }
};

} // namespace ligero::net

