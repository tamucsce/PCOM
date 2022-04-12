#include <chrono>

#include <boost/program_options.hpp>

#include <ligero/tcp/server.hpp>
#include <ligero/tcp/client.hpp>
#include <ligero/network/serialize.hpp>

using namespace ligero::net;
using namespace std::chrono_literals;
using namespace std::string_literals;

using serializer_t = serializer<boost_portable_binary>;
using tcp_server = basic_tcp_server<serializer_t>;
using tcp_client = basic_tcp_client<serializer_t>;

void run_server(std::chrono::seconds accept_time) {

    /*  Initialize a executor */
    asio::io_context executor;

    /*  Construct a TCP server */
    tcp_server sock(executor);

    /*  Optional: set worker threads */
    sock.threads(2);

    DEBUG << "Binding to 127.0.0.1:5555";

    /*  Bind the server to localhost:5555 */
    auto accepted = sock.bind("127.0.0.1", 5555, accept_time);

    INFO << "Accepted " << accepted << " connections";

    /*  Notify clients to start */
    sock.broadcast("start"s);

    /*  Deserialize received messages and send back */
    {
        auto [msgs, kickedout] = sock.collect(1s);

        if (kickedout) {
            WARNING << kickedout << " clients have been kicked out";
        }
        
        int sum = 0;
        for (auto& msg : msgs) {
            /*  Deserialize using server's archive. Alternatively, one could write:
             *  
             *  int a, b, c;
             *  sock.archive().unpack(msg, a, b, c);
             */
            auto [n] = sock.archive().unpack<int>(msg);
            sum += n;
        }
        DEBUG << "Sum is: " << sum;
        sock.broadcast(sum);
    }

    INFO << "Done!";
}

// ------------------------------------------------------------

void run_client(int input) {

    /*  Additional thread hint, optimized for single thread */
    asio::io_context executor(1);

    tcp_client sock(executor);

    /*  Connect to remote server  */
    sock.connect("127.0.0.1", "5555");

    /*  Wait until server is ready... */
    auto [str] = sock.receive<std::string>(10s);

    DEBUG << "Server says: " << str;

    {
        sock.send(input);
        auto [sum] = sock.receive<int>(10s);

        INFO << "Get sum: " << sum;
    }

    auto [up, down] = sock.communication_cost(storage_unit::Byte);

    INFO << "Uploaded: " << up << " Bytes, "
         << "Downloaded: " << down << " Bytes";

    INFO << "Done!";
}

// ------------------------------------------------------------

int main(int argc, char *argv[]) {
    namespace po = boost::program_options;

    po::options_description desc("Options");
    desc.add_options()
        ("client,c",   "Act like client")
        ("server,s",   "Act like server")
        ("duration,d", po::value<int>(), "Registration time")
        ("input,i",    po::value<int>(), "Input value");

    try {
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);

        int duration_sec = 5, input = 1;

        if (vm.count("client")) {
            if (vm.count("input")) { input = vm["input"].as<int>(); }
            run_client(input);
        }
        else if (vm.count("server")) {
            if (vm.count("duration")) { duration_sec = vm["duration"].as<int>(); }
            run_server(std::chrono::seconds(duration_sec));
        }
        else {
            std::cout << desc << std::endl;
        }
    }
    catch (...) {
        ERROR << boost::current_exception_diagnostic_information();
        return 1;
    }

    return 0;
}
