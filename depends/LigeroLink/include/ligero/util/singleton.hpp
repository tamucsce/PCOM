#pragma once

namespace ligero {

template <class T>
struct singleton : private T {
    static T& instance() {
        static singleton<T> s;  // compiler with c++11 will ensure thread safety
        return s;
    }
};

}  // namespace ligero

