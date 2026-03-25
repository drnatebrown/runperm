#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "orbit/common.hpp"
#include <cmath>
#include <cassert>
#include <optional>
#include <iostream>

namespace orbit {

inline constexpr std::optional<double> DEFAULT_LENGTH_CAPPING = 8.0;
inline constexpr std::optional<ulint> DEFAULT_BALANCING = 16;

struct split_params {
    std::optional<double> length_capping;
    std::optional<ulint> balancing;

    split_params() : length_capping(DEFAULT_LENGTH_CAPPING), balancing(DEFAULT_BALANCING) {}
    split_params(std::optional<double> length_capping, std::optional<ulint> balancing)
    : length_capping(std::move(length_capping)), balancing(std::move(balancing)) {}

    bool operator==(const split_params& other) const {
        return length_capping == other.length_capping && balancing == other.balancing;
    }
    bool operator!=(const split_params& other) const { return !(*this == other); }

    size_t serialize(std::ostream& out) const {
        size_t written_bytes = 0;
        bool has_length_capping = length_capping.has_value();
        out.write(reinterpret_cast<const char*>(&has_length_capping), sizeof(has_length_capping));
        written_bytes += sizeof(has_length_capping);
        if (length_capping.has_value()) {
            double length_capping_value = length_capping.value();
            out.write(reinterpret_cast<const char*>(&length_capping_value), sizeof(length_capping_value));
            written_bytes += sizeof(length_capping_value);
        }
        bool has_balancing = balancing.has_value();
        out.write(reinterpret_cast<const char*>(&has_balancing), sizeof(has_balancing));
        written_bytes += sizeof(has_balancing);
        if (balancing.has_value()) {
            ulint balancing_value = balancing.value();
            out.write(reinterpret_cast<const char*>(&balancing_value), sizeof(balancing_value));
            written_bytes += sizeof(balancing_value);
        }
        return written_bytes;
    }

    void load(std::istream& in) {
        length_capping = std::nullopt;
        balancing = std::nullopt;
        bool has_length_capping;
        in.read(reinterpret_cast<char*>(&has_length_capping), sizeof(has_length_capping));
        if (has_length_capping) {
            double length_capping_value;
            in.read(reinterpret_cast<char*>(&length_capping_value), sizeof(length_capping_value));
            length_capping = length_capping_value;
        }
        bool has_balancing;
        in.read(reinterpret_cast<char*>(&has_balancing), sizeof(has_balancing));
        if (has_balancing) {
            ulint balancing_value;
            in.read(reinterpret_cast<char*>(&balancing_value), sizeof(balancing_value));
            balancing = balancing_value;
        }
    }
};

inline split_params NO_SPLITTING = split_params(std::nullopt, std::nullopt);
inline split_params DEFAULT_SPLITTING = split_params(DEFAULT_LENGTH_CAPPING, DEFAULT_BALANCING);
inline split_params ONLY_LENGTH_CAPPING = split_params(DEFAULT_LENGTH_CAPPING, std::nullopt);
inline split_params ONLY_BALANCING = split_params(std::nullopt, DEFAULT_BALANCING);

template<class int_vector_t>
struct split_result {
    int_vector_t lengths;
    int_vector_t img_rank_inv;
    ulint max_length;
};

template<class int_vector_t>
inline void split_by_length_capping(
    const int_vector_t& lengths, 
    const int_vector_t& img_rank_inv, 
    const ulint domain, 
    const double length_capping_factor, 
    split_result<int_vector_t>& result
) {
    assert(lengths.size() == img_rank_inv.size());
    assert(length_capping_factor > 0.0);

    double avg_run_length = static_cast<double>(domain) / static_cast<double>(lengths.size());
    ulint desired_max_allowed_length = static_cast<ulint>(std::ceil(avg_run_length * length_capping_factor));
    uchar bits = bit_width(desired_max_allowed_length);
    ulint max_allowed_length = max_val(bits);

    size_t new_intervals_upper_bound = std::ceil(static_cast<double>(lengths.size()) / length_capping_factor);

    split_by_max_allowed_length(lengths, img_rank_inv, domain, max_allowed_length, result, new_intervals_upper_bound);
}

template<class int_vector_t>
inline void split_by_max_allowed_length(
    const int_vector_t& lengths, 
    const int_vector_t& img_rank_inv,
    const ulint domain,
    const ulint max_allowed_length, 
    split_result<int_vector_t>& result,
    const size_t new_intervals_upper_bound = 0
) {
    assert(lengths.size() == img_rank_inv.size());
    assert(max_allowed_length > 0);
    size_t intervals_after_splitting_upper_bound = 
        (new_intervals_upper_bound == 0) ? domain - 1 : lengths.size() + new_intervals_upper_bound - 1;

    // For each original run, how many new intervals were added up to but not including this run?
    int_vector_t input_splits_exclusive_cumsum(lengths.size(), bit_width(intervals_after_splitting_upper_bound));
    
    // First pass to determine the number of intervals after splitting in input order, 
    // the number of new intervals added up to but not including each run, 
    // and the max length of the new intervals
    size_t num_intervals_after_splitting = 0;
    size_t cumulative_new_intervals = 0;
    result.max_length = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        input_splits_exclusive_cumsum[i] = cumulative_new_intervals;

        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                remaining -= chunk;
                ++num_intervals_after_splitting;
                ++cumulative_new_intervals;
            }
            --cumulative_new_intervals;
            result.max_length = max_allowed_length;
        } else {
            result.max_length = std::max(result.max_length, lengths[i]);
            ++num_intervals_after_splitting;
        }
    }

    result.lengths = int_vector_t(num_intervals_after_splitting, bit_width(result.max_length));
    result.img_rank_inv = int_vector_t(num_intervals_after_splitting, bit_width(num_intervals_after_splitting - 1));

    // Second pass to fill the lengths and img_rank_inv arrays, in output order
    size_t curr_img_rank_inv_idx = 0;
    for (size_t i = 0; i < lengths.size(); ++i) {
        size_t j = img_rank_inv[i];
        size_t length = lengths[j];
        size_t num_splits = 0;
        if (length > max_allowed_length) {
            ulint remaining = length;
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                remaining -= chunk;
                result.lengths[j + input_splits_exclusive_cumsum[j] + num_splits] = chunk;
                ++num_splits;
            }
            --num_splits;
        } else {
            result.lengths[j + input_splits_exclusive_cumsum[j]] = length;
        }
        
        // Fill the img_rank_inv array
        size_t curr_img_rank_inv_val = j + input_splits_exclusive_cumsum[j];
        for (size_t k = 0; k < num_splits + 1; ++k) {
            result.img_rank_inv[curr_img_rank_inv_idx] = curr_img_rank_inv_val + k;
            ++curr_img_rank_inv_idx;
        }
    }
}

template<class int_vector_t>
inline void split_by_balancing(const int_vector_t& lengths, const int_vector_t& img_rank_inv, const ulint domain, const ulint balancing_factor, split_result<int_vector_t>& result) {
    assert(lengths.size() == img_rank_inv.size());
    // TODO
    result.lengths = lengths;
    result.img_rank_inv = img_rank_inv;
    // result.max_length = 0;
}

} // namespace orbit

#endif