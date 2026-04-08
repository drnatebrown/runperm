#ifndef _INTERVAL_ENCODING_IMPL_HPP
#define _INTERVAL_ENCODING_IMPL_HPP

#include "orbit/common.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include "orbit/internal/ds/packed_vector_aligned.hpp"
#include "orbit/internal/move/move_splitting.hpp"

#include <optional>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace orbit {

/* PERMUTATION UTILITIES */
template<typename container_t>
inline std::vector<ulint> get_inverse_permutation(const container_t& permutation) {
    std::vector<ulint> inverse_permutation(permutation.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        inverse_permutation[permutation[i]] = i;
    }
    return inverse_permutation;
}

inline std::pair<std::vector<ulint>, std::vector<ulint>> get_permutation_intervals(const std::vector<ulint> &permutation, ulint* max_length_ret = nullptr) {
    std::vector<ulint> lengths;
    std::vector<ulint> images;
    ulint max_length = 0;
    for (size_t i = 0; i < permutation.size(); ++i) {
        if (i == 0 || permutation[i] != permutation[i - 1] + 1) {
            if (!lengths.empty()) {
                max_length = std::max(max_length, lengths.back());
            }
            lengths.push_back(1);
            images.push_back(permutation[i]);
        } else {
            ++lengths.back();
        }
    }
    if (!lengths.empty()) {
        max_length = std::max(max_length, lengths.back());
    }
    if (max_length_ret) {
        *max_length_ret = max_length;
    }
    return {lengths, images};
}

template<typename container_t, typename int_vector_t = int_vector_aligned>
inline std::tuple<int_vector_t, ulint> starts_to_lengths(const container_t& starts, const ulint domain) {
    int_vector_t lengths(starts.size(), bit_width(domain - 1));
    ulint max_length = 0;
    for (size_t i = 0; i < lengths.size() - 1; ++i) {
        lengths[i] = starts[i + 1] - starts[i];
        max_length = std::max(max_length, static_cast<ulint>(lengths[i]));
    }
    lengths[lengths.size() - 1] = domain - starts[starts.size() - 1];
    max_length = std::max(max_length, static_cast<ulint>(lengths[lengths.size() - 1]));
    return {lengths, max_length};
}

template<typename container_t>
inline std::tuple<ulint, ulint> sum_and_max(const container_t& data) {
    ulint sum = 0;
    ulint max = 0;
    for (const auto& element : data) {
        sum += element;
        max = std::max(max, static_cast<ulint>(element));
    }
    return {sum, max};
}

template<typename int_vector_t = int_vector_aligned, typename container_t>
inline int_vector_t compute_img_rank_inv(const container_t& interval_output_starts) {
    int_vector_t img_rank_inv(interval_output_starts.size(), bit_width(interval_output_starts.size() - 1));
    std::iota(img_rank_inv.begin(), img_rank_inv.end(), 0);
    std::sort(img_rank_inv.begin(), img_rank_inv.end(), [&](ulint a, ulint b) { 
        return static_cast<ulint>(interval_output_starts[a]) < static_cast<ulint>(interval_output_starts[b]); 
    });
    return img_rank_inv;
}

template<typename int_vector_t = int_vector_aligned>
class interval_encoding_impl {
public:

    interval_encoding_impl() = default;
    interval_encoding_impl(const std::vector<ulint>& permutation, const split_params& sp = split_params{}) {
        *this = from_permutation(permutation, sp);
    }
    interval_encoding_impl(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const split_params& sp = split_params{}) {
        *this = from_lengths_and_images(lengths, images, sp);
    }
    interval_encoding_impl(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const ulint domain, const ulint max_length, const split_params& sp = split_params{}) {
        *this = from_lengths_and_images(lengths, images, domain, max_length, sp);
    }

    static interval_encoding_impl<int_vector_t> from_permutation(const std::vector<ulint>& permutation, const split_params& sp = split_params{}) {
        ulint max_length = 0;
        auto [lengths, images] = get_permutation_intervals(permutation, &max_length);
        assert(lengths.size() == images.size());

        return from_lengths_and_images(lengths, images, permutation.size(), max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_img_rank_inv(const container1_t& lengths, const container2_t& img_rank_inv, const split_params& sp = split_params{}) {
        assert(lengths.size() == img_rank_inv.size());
        auto [domain, max_length] = sum_and_max(lengths);
        
        return from_lengths_and_img_rank_inv(lengths, img_rank_inv, domain, max_length, sp);
    }
    
    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_img_rank_inv(const container1_t& lengths, const container2_t& img_rank_inv, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(lengths.size() == img_rank_inv.size());
        interval_encoding_impl<int_vector_t> enc;
        enc.set_initial_values(domain, lengths.size(), max_length, sp);
        enc.init_img_rank_inv(lengths, img_rank_inv);
        return enc;
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_img_rank(const container1_t& lengths, const container2_t& img_rank, const split_params& sp = split_params{}) {
        assert(lengths.size() == img_rank.size());

        auto [domain, max_length] = sum_and_max(lengths);
        return from_lengths_and_img_rank(lengths, img_rank, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_img_rank(const container1_t& lengths, const container2_t& img_rank, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(lengths.size() == img_rank.size());

        interval_encoding_impl<int_vector_t> enc;
        enc.set_initial_values(domain, lengths.size(), max_length, sp);
        enc.init_img_rank(lengths, img_rank);
        return enc;
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_images(const container1_t& lengths, const container2_t& images, const split_params& sp = split_params{}) {
        assert(lengths.size() == images.size());

        auto [domain, max_length] = sum_and_max(lengths);
        return from_lengths_and_images(lengths, images, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_lengths_and_images(const container1_t& lengths, const container2_t& images, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(lengths.size() == images.size());

        int_vector_t img_rank_inv = compute_img_rank_inv<int_vector_t>(images);
        return from_lengths_and_img_rank_inv(lengths, img_rank_inv, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_img_rank_inv(const container1_t& starts, const container2_t& img_rank_inv, const size_t domain, const split_params& sp = split_params{}) {
        assert(starts.size() == img_rank_inv.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_img_rank_inv(lengths, img_rank_inv, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_img_rank_inv(const container1_t& starts, const container2_t& img_rank_inv, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(starts.size() == img_rank_inv.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_img_rank_inv(lengths, img_rank_inv, domain, calculated_max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_img_rank(const container1_t& starts, const container2_t& img_rank, const size_t domain, const split_params& sp = split_params{}) {
        assert(starts.size() == img_rank.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_img_rank(lengths, img_rank, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_img_rank(const container1_t& starts, const container2_t& img_rank, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(starts.size() == img_rank.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_img_rank(lengths, img_rank, domain, calculated_max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_images(const container1_t& starts, const container2_t& images, const size_t domain, const split_params& sp = split_params{}) {
        assert(starts.size() == images.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_images(lengths, images, domain, max_length, sp);
    }

    template<typename container1_t, typename container2_t>
    static interval_encoding_impl<int_vector_t> from_starts_and_images(const container1_t& starts, const container2_t& images, const size_t domain, const ulint max_length, const split_params& sp = split_params{}) {
        assert(starts.size() == images.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_images(lengths, images, domain, calculated_max_length, sp);
    }

    template<typename T, typename container_t>
    std::vector<T> split_run_data_with_copy(const container_t& original_lengths, const std::vector<T>& run_data) {
        assert(original_lengths.size() == run_data.size());
        assert(original_lengths.size() == this->runs_);

        std::vector<T> new_run_data(this->intervals_);
        size_t new_interval_idx = 0;
        for (size_t i = 0; i < original_lengths.size(); ++i) {
            ulint remaining = original_lengths[i];
            while (remaining > 0) {
                ulint chunk = get_length(new_interval_idx);
                new_run_data[new_interval_idx] = run_data[i];
                remaining -= chunk;
                ++new_interval_idx;
            }
        }
        assert(new_interval_idx == this->intervals_);

        return new_run_data;
    }

    ulint get_length(size_t i) const {
        return lengths[i];
    }
    ulint get_img_rank_inv(size_t i) const {
        return img_rank_inv[i];
    }
    size_t domain() const {
        return domain_;
    }
    size_t runs() const {
        return runs_;
    }
    size_t intervals() const {
        return intervals_;
    }
    ulint max_length() const {
        return max_length_;
    }
    split_params get_split_params() const {
        return split_params_;
    }

protected:
    int_vector_t lengths;
    int_vector_t img_rank_inv;
    split_params split_params_;
    size_t domain_;
    size_t runs_;
    size_t intervals_;
    ulint max_length_;

    void set_initial_values(size_t domain, size_t runs, ulint max_length, const split_params& split_params) {
        this->domain_ = domain;
        this->runs_ = runs;
        this->max_length_ = max_length;
        this->split_params_ = split_params;
    }
    
    template<typename container1_t, typename container2_t>
    void init_img_rank(const container1_t& lengths, const container2_t& img_rank) {
        uchar length_bits = bit_width(this->max_length_);
        uchar img_rank_inv_bits = bit_width(this->runs_ - 1);

        int_vector_t curr_lengths(this->runs_, length_bits);
        int_vector_t curr_img_rank_inv(this->runs_, img_rank_inv_bits);
        for (size_t i = 0; i < this->runs_; ++i) {
            curr_lengths[i] = lengths[i];
            curr_img_rank_inv[img_rank[i]] = i;
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_img_rank_inv, new_max_length);

        this->lengths = std::move(curr_lengths);
        this->img_rank_inv = std::move(curr_img_rank_inv);
        this->max_length_ = new_max_length;
        this->intervals_ = this->lengths.size();
    }

    template<typename container1_t, typename container2_t>
    void init_img_rank_inv(const container1_t& lengths, const container2_t& img_rank_inv) { 
        uchar length_bits = bit_width(this->max_length_);
        uchar img_rank_inv_bits = bit_width(this->runs_ - 1);

        int_vector_t curr_lengths(this->runs_, length_bits);
        int_vector_t curr_img_rank_inv(this->runs_, img_rank_inv_bits);
        for (size_t i = 0; i < this->runs_; ++i) {
            curr_lengths[i] = lengths[i];
            curr_img_rank_inv[i] = img_rank_inv[i];
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_img_rank_inv, new_max_length);

        this->lengths = std::move(curr_lengths);
        this->img_rank_inv = std::move(curr_img_rank_inv);
        this->max_length_ = new_max_length;
        this->intervals_ = this->lengths.size();
    }

    void apply_splitting(int_vector_t& curr_lengths, int_vector_t& curr_img_rank_inv, ulint& new_max_length) {
        if (this->split_params_ == NO_SPLITTING) {
            new_max_length = this->max_length_;
            return;
        }
        
        split_result<int_vector_t> split_result;
        
        if (this->split_params_.length_capping) {
            split_by_length_capping(curr_lengths, curr_img_rank_inv, this->domain_, *this->split_params_.length_capping, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_img_rank_inv = std::move(split_result.img_rank_inv);
            new_max_length = split_result.max_length;
        }
        if (this->split_params_.balancing) {
            split_by_balancing(curr_lengths, curr_img_rank_inv, this->domain_, *this->split_params_.balancing, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_img_rank_inv = std::move(split_result.img_rank_inv);
            new_max_length = split_result.max_length;
        }
    }
};

} // namespace orbit

#endif // end include guard _INTERVAL_ENCODING_IMPL_HPP