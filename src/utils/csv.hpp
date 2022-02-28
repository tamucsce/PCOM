#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>

namespace psi::util {

namespace fs = std::filesystem;

std::fstream open_file(const fs::path& path) {
    if (!fs::exists(path)) {
        throw std::runtime_error("File does not exist: " + path.string());
    }

    std::fstream fs(path);
    if (!fs.is_open()) {
        throw std::runtime_error("Cannot open file: " + path.string());
    }

    return fs;
}

template <typename T>
std::vector<T> read_line(std::fstream& is, char delimiter = '\n') {
    std::string str;
    std::vector<T> result;
        
    while (std::getline(is, str, delimiter)) {
	if (str.empty())
	    continue;
        std::stringstream ss(str);
        T value;
        ss >> value;
        result.push_back(std::move(value));
    }
    return result;
}

template <typename Iter>
std::ostream& write_line(std::fstream& os,
                         Iter begin, Iter end,
                         char delimiter = '\n')
{
    assert(std::distance(begin, end) > 0);
    for (auto it = begin; it != std::prev(end); it++) {
        os << *it << delimiter;
    }
    os << *std::prev(end);
    os << std::endl;
    return os;
}

}
