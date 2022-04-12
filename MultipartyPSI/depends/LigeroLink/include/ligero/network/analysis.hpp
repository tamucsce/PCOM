#pragma once

#include <numeric>
#include <tuple>
#include <vector>

namespace ligero::net {

enum class storage_unit : unsigned short {
    Byte = 0,
    KB,
    MB,
    GB,
    TB,
};

struct communication_analysis {

    void analysis_emplace_sent(size_t size) { sent_.emplace_back(size); }
    void analysis_emplace_received(size_t size) { received_.emplace_back(size); }

    const std::vector<size_t>& uploaded_bytes() const { return sent_; }
    const std::vector<size_t>& downloaded_bytes() const { return received_; }

    std::pair<size_t, size_t> communication_cost(storage_unit unit = storage_unit::KB) const {
        size_t upload = std::reduce(sent_.cbegin(), sent_.cend());
        size_t download = std::reduce(received_.cbegin(), received_.cend());

        using Unit = storage_unit;
        size_t measure = 1;
        constexpr size_t kb = 1024;
        
        switch (unit) {
        case Unit::Byte:
            break;
        case Unit::KB:
            measure = kb; break;
        case Unit::MB:
            measure = kb * kb; break;
        case Unit::GB:
            measure = kb * kb * kb; break;
        case Unit::TB:
            measure = kb * kb * kb * kb; break;
        default:
            throw std::runtime_error("Unexpected measurement");
        }
        return { upload / measure, download / measure };
    }

protected:
    std::vector<size_t> sent_, received_;
};

}
