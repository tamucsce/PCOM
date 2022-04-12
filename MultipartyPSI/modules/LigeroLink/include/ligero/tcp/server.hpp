#pragma once

#include <algorithm>
#include <chrono>
#include <list>
#include <span>

#include <ligero/exception.hpp>
#include <ligero/util/log.hpp>
#include <ligero/util/timer.hpp>

#include <boost/asio.hpp>
#include <boost/beast/core.hpp>

namespace ligero::net {

using namespace boost;
using namespace std::chrono_literals;

using asio::awaitable;
using asio::co_spawn;
using asio::detached;
using asio::use_awaitable;
using asio::ip::tcp;

constexpr auto version = "aya-websocket-0.3";

template <typename SerializePolicy>
struct basic_tcp_server {

    using socket_t = beast::tcp_stream;

protected:
    size_t num_threads_ = std::thread::hardware_concurrency();
    SerializePolicy archive_;
    asio::io_context& executor_;
    std::list<socket_t> sockets_;
    std::list<beast::flat_buffer> buffers_;

public:

    basic_tcp_server(asio::io_context& exe) : executor_(exe) { }
    basic_tcp_server(asio::io_context& exe, size_t allocate_hint)
        : archive_(allocate_hint), executor_(exe) { }

    // Disable copy construction
    basic_tcp_server(const basic_tcp_server&) = delete;
    basic_tcp_server& operator=(const basic_tcp_server&) = delete;

    virtual ~basic_tcp_server() = default;

    size_t threads() const noexcept { return num_threads_; }
    void threads(size_t t) noexcept { num_threads_ = t; }

    auto& archive() { return archive_; }
    const auto& archive() const { return archive_; }
    auto& executor() noexcept { return executor_; }
    auto& sockets() noexcept { return sockets_; }
    auto& buffers() noexcept { return buffers_; }
    
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

    /** Initialize the communication channel and accept incoming connection.
     *  Must be called before performing any other operations.
     */
    size_t bind(std::string endpoint, uint16_t port, std::chrono::nanoseconds bind_timeout) {
        tcp::acceptor acceptor{ executor_, { asio::ip::make_address(endpoint), port }};

        /*  Spawn a timeout timer to cancel TCP accept */
        asio::high_resolution_timer timer(executor_);
        timer.expires_after(bind_timeout);
        
        co_spawn(executor_,
                 [&timer, &acceptor]() -> awaitable<void> {
                     co_await timer.async_wait(use_awaitable);
                     acceptor.cancel();
                 }, detached);

        /*  Spawn the accept coroutine */
        co_spawn(executor_,
                 [&acceptor, this]() -> awaitable<void> {
                     auto executor = co_await asio::this_coro::executor;
                     auto strand = asio::make_strand(executor);
                     
                     for (;;) {
                         tcp::socket socket = co_await coro_tcp_accept(acceptor);
                         sockets_.emplace_back(std::move(socket));
                     }
                 }, detached);

        /*  Block wait the coroutine finish */
        launch_coroutine();

        buffers_.resize(sockets_.size());
        
        return sockets_.size();
    }

    /** Serialize and broadcast a (sequence of) message to connected peers
     *  in parallel.
     */
    template <typename... Args>
    size_t broadcast(Args&&... data) {
        auto __t = util::make_timer("Ligero", "Net", "Broadcast");

        auto view = archive_.pack(std::forward<Args>(data)...);
        auto buf = asio::buffer(reinterpret_cast<const void*>(view.data()), view.size());

        uint64_t message_size = view.size();
        auto header_buffer = asio::buffer(reinterpret_cast<const char*>(&message_size), sizeof(uint64_t));

        for (auto& socket : sockets_) {
            co_spawn(executor_,
                     [&]() -> awaitable<void> {
                         co_await asio::async_write(socket, header_buffer, use_awaitable);
                         co_await asio::async_write(socket, buf, use_awaitable);
                     },
                     detached);
        }

        launch_coroutine();
        return view.size();
    }

    std::pair<std::vector<std::span<const char>>, size_t>
    collect(std::chrono::nanoseconds timeout) {
        auto __t = util::make_timer("Ligero", "Net", "Collect");

        std::vector<std::exception_ptr> exceptions(sockets_.size());
        size_t kicked_out = 0;

        /*  Set receive timeout */
        expire_after(timeout);

        auto buffer_it = buffers_.begin();
        for (auto it = sockets_.begin(); it != sockets_.end(); it++, buffer_it++) {
            auto pos = std::distance(sockets_.begin(), it);

            co_spawn(executor_,
                     coro_receive(*it, *buffer_it),
                     [&, pos](std::exception_ptr e) { exceptions[pos] = e; });
        }
        launch_coroutine();

        /*  Remember to reset timeout */
        reset_timeout();

        /*  Mark bad connections, if any */
        std::vector<decltype(sockets_.begin())> socket_iters;
        std::vector<decltype(buffers_.begin())> buffer_iters;
        for (auto i = 0; i < exceptions.size(); i++) {
            if (exceptions[i]) {
                try {
                    std::rethrow_exception(exceptions[i]);
                }
                catch (const system::system_error& e) {
                    WARNING << "Socket [" << i <<  "] failed during collect: " << e.what();
                    socket_iters.emplace_back(std::next(sockets_.begin(), i));
                    buffer_iters.emplace_back(std::next(buffers_.begin(), i));
                    ++kicked_out;
                }
            }
        }

        /*  Close the bad connection and erase from socket list */
        std::for_each(socket_iters.begin(), socket_iters.end(),
                      [this](auto&& it) { sockets_.erase(it); });
        std::for_each(buffer_iters.begin(), buffer_iters.end(),
                      [this](auto&& it) { buffers_.erase(it); });
        assert(buffers_.size() == sockets_.size());

        if (sockets_.empty()) {
            LIGERO_THROW( network_error()
                          << throw_reason("No message received")
                          << throw_error_code(error_code::no_message_received) );
        }

        /* Encapsulate return type to make life easier */
        std::vector<std::span<const char>> views;
        for (auto& buf : buffers_) {
            views.emplace_back(reinterpret_cast<const char*>(buf.data().data()),
                               buf.data().size());
        }
        
        return std::make_pair(views, kicked_out);
    }

    // ------------------------------------------------------------
    
protected:
    /** A tcp connection can be established after two handshakes:
     *    1. Accept TCP handshake
     *    2. Accept Websocket handshake
     *  We assume a client will perform websocket handshake right after TCP
     *  connection (i.e. no significant delay).
     **/

    /** Set a high resolution timer */
    void expire_after(std::chrono::nanoseconds duration) {
        for (auto& socket : sockets_) {
            beast::get_lowest_layer(socket).expires_after(duration);
        }
    }

    /** Reset all timers */
    void reset_timeout() {
        for (auto& socket : sockets_) {
            beast::get_lowest_layer(socket).expires_never();
        }
    }

    awaitable<tcp::socket> coro_tcp_accept(tcp::acceptor& acceptor) {
        auto executor = co_await asio::this_coro::executor;
        
        tcp::socket socket = co_await acceptor.async_accept(executor, use_awaitable);
        socket.set_option(asio::socket_base::keep_alive{true});
        socket.set_option(tcp::no_delay{true});
        
        co_return std::move(socket);
    }

    template <typename Buffer>
    awaitable<void> async_read(socket_t& socket, Buffer& buffer, size_t bytes) {
        auto buf = buffer.prepare(bytes);
        auto n = co_await asio::async_read(socket, buf, use_awaitable);
        buffer.commit(n);
    }

    awaitable<void> coro_receive(socket_t& socket, beast::flat_buffer& buffer) {
        buffer.consume(buffer.size());
        
        constexpr size_t header_size = 8;
        uint64_t message_len = 0;
        
        co_await async_read(socket, buffer, header_size);

        std::copy(reinterpret_cast<const char*>(buffer.data().data()),
                  reinterpret_cast<const char*>(buffer.data().data()) + buffer.data().size(),
                  reinterpret_cast<char*>(&message_len));

        buffer.consume(buffer.size());

        co_await async_read(socket, buffer, message_len);
    }
};

} // namespace ligero::net
