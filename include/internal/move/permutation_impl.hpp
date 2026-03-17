#ifndef _PERMUTATION_HPP
#define _PERMUTATION_HPP

#include "common.hpp"
#include "internal/ds/packed_vector.hpp"
#include "internal/ds/packed_vector_aligned.hpp"
#include "internal/move/move_splitting.hpp"

#include <optional>
#include <numeric>
#include <algorithm>
#include <cassert>

/* PERMUTATION UTILITIES */
template<class Container>
inline std::vector<ulint> get_inverse_permutation(const Container& permutation) {
    std::vector<ulint> inverse_permutation(permutation.size());
    for (size_t i = 0; i < permutation.size(); ++i) {
        inverse_permutation[permutation[i]] = i;
    }
    return inverse_permutation;
}

inline std::pair<std::vector<ulint>, std::vector<ulint>> get_permutation_intervals(const std::vector<ulint> &permutation, ulint* max_length_ret = nullptr) {
    std::vector<ulint> lengths;
    std::vector<ulint> interval_permutation;
    ulint max_length = 0;
    for (size_t i = 0; i < permutation.size(); ++i) {
        if (i == 0 || permutation[i] != permutation[i - 1] + 1) {
            if (!lengths.empty()) {
                max_length = std::max(max_length, lengths.back());
            }
            lengths.push_back(1);
            interval_permutation.push_back(permutation[i]);
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
    return {lengths, interval_permutation};
}

template<class Container, class IntVectorType = IntVectorAligned>
inline std::tuple<IntVectorType, ulint> starts_to_lengths(const Container& starts, const ulint domain) {
    IntVectorType lengths(starts.size(), bit_width(domain - 1));
    ulint max_length = 0;
    for (size_t i = 0; i < lengths.size() - 1; ++i) {
        lengths[i] = starts[i + 1] - starts[i];
        max_length = std::max(max_length, static_cast<ulint>(lengths[i]));
    }
    lengths[lengths.size() - 1] = domain - starts[starts.size() - 1];
    max_length = std::max(max_length, static_cast<ulint>(lengths[lengths.size() - 1]));
    return {lengths, max_length};
}

template<class Container>
inline std::tuple<ulint, ulint> sum_and_max(const Container& data) {
    ulint sum = 0;
    ulint max = 0;
    for (const auto& element : data) {
        sum += element;
        max = std::max(max, static_cast<ulint>(element));
    }
    return {sum, max};
}

template<class IntVectorType = IntVectorAligned, class Container>
inline IntVectorType compute_tau_inv(const Container& interval_output_starts) {
    IntVectorType tau_inv(interval_output_starts.size(), bit_width(interval_output_starts.size() - 1));
    std::iota(tau_inv.begin(), tau_inv.end(), 0);
    std::sort(tau_inv.begin(), tau_inv.end(), [&](ulint a, ulint b) { 
        return static_cast<ulint>(interval_output_starts[a]) < static_cast<ulint>(interval_output_starts[b]); 
    });
    return tau_inv;
}

template<class IntVectorType = IntVectorAligned>
class PermutationImpl {
public:

    PermutationImpl() = default;

    PermutationImpl(const std::vector<ulint>& permutation, const SplitParams& split_params = SplitParams()) {
        *this = from_permutation(permutation, split_params);
    }
    PermutationImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const SplitParams& split_params = SplitParams()) {
        *this = from_lengths_and_interval_permutation(lengths, interval_permutation, split_params);
    }
    PermutationImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        *this = from_lengths_and_interval_permutation(lengths, interval_permutation, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_permutation(const std::vector<ulint>& permutation, const SplitParams& split_params = SplitParams()) {
        ulint max_length = 0;
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation, &max_length);
        assert(lengths.size() == interval_permutation.size());

        return from_lengths_and_interval_permutation(lengths, interval_permutation, permutation.size(), max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_tau_inv(const Container1& lengths, const Container2& tau_inv, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau_inv.size());
        auto [domain, max_length] = sum_and_max(lengths);
        
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }
    
    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_tau_inv(const Container1& lengths, const Container2& tau_inv, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau_inv.size());
        PermutationImpl<IntVectorType> permutation;
        permutation.set_initial_values(domain, lengths.size(), max_length, split_params);
        permutation.init_tau_inv(lengths, tau_inv);
        return permutation;
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_tau(const Container1& lengths, const Container2& tau, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau.size());

        auto [domain, max_length] = sum_and_max(lengths);
        return from_lengths_and_tau(lengths, tau, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_tau(const Container1& lengths, const Container2& tau, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau.size());

        PermutationImpl<IntVectorType> permutation;
        permutation.set_initial_values(domain, lengths.size(), max_length, split_params);
        permutation.init_tau(lengths, tau);
        return permutation;
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_interval_permutation(const Container1& lengths, const Container2& interval_permutation, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        auto [domain, max_length] = sum_and_max(lengths);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_lengths_and_interval_permutation(const Container1& lengths, const Container2& interval_permutation, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        IntVectorType tau_inv = compute_tau_inv<IntVectorType>(interval_permutation);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_tau_inv(const Container1& starts, const Container2& tau_inv, const size_t domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau_inv.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_tau_inv(const Container1& starts, const Container2& tau_inv, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau_inv.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, calculated_max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_tau(const Container1& starts, const Container2& tau, const size_t domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_tau(lengths, tau, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_tau(const Container1& starts, const Container2& tau, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_tau(lengths, tau, domain, calculated_max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_interval_permutation(const Container1& starts, const Container2& interval_permutation, const size_t domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == interval_permutation.size());
        auto [lengths, max_length] = starts_to_lengths(starts, domain);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, max_length, split_params);
    }

    template<class Container1, class Container2>
    static PermutationImpl<IntVectorType> from_starts_and_interval_permutation(const Container1& starts, const Container2& interval_permutation, const size_t domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == interval_permutation.size());
        auto [lengths, calculated_max_length] = starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, calculated_max_length, split_params);
    }

    template<typename T, class Container>
    std::vector<T> split_run_data_with_copy(const Container& original_lengths, const std::vector<T>& run_data) {
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
    ulint get_tau_inv(size_t i) const {
        return tau_inv[i];
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
    SplitParams get_split_params() const {
        return split_params_;
    }

protected:
    IntVectorType lengths;
    IntVectorType tau_inv;
    SplitParams split_params_;
    size_t domain_;
    size_t runs_;
    size_t intervals_;
    ulint max_length_;

    void set_initial_values(size_t domain, size_t runs, ulint max_length, const SplitParams& split_params) {
        this->domain_ = domain;
        this->runs_ = runs;
        this->max_length_ = max_length;
        this->split_params_ = split_params;
    }
    
    template<class Container1, class Container2>
    void init_tau(const Container1& lengths, const Container2& tau) {
        uchar length_bits = bit_width(this->max_length_);
        uchar tau_inv_bits = bit_width(this->runs_ - 1);

        IntVectorType curr_lengths(this->runs_, length_bits);
        IntVectorType curr_tau_inv(this->runs_, tau_inv_bits);
        for (size_t i = 0; i < this->runs_; ++i) {
            curr_lengths[i] = lengths[i];
            curr_tau_inv[tau[i]] = i;
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_tau_inv, new_max_length);

        this->lengths = std::move(curr_lengths);
        this->tau_inv = std::move(curr_tau_inv);
        this->max_length_ = new_max_length;
        this->intervals_ = this->lengths.size();
    }

    template<class Container1, class Container2>
    void init_tau_inv(const Container1& lengths, const Container2& tau_inv) { 
        uchar length_bits = bit_width(this->max_length_);
        uchar tau_inv_bits = bit_width(this->runs_ - 1);

        IntVectorType curr_lengths(this->runs_, length_bits);
        IntVectorType curr_tau_inv(this->runs_, tau_inv_bits);
        for (size_t i = 0; i < this->runs_; ++i) {
            curr_lengths[i] = lengths[i];
            curr_tau_inv[i] = tau_inv[i];
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_tau_inv, new_max_length);

        this->lengths = std::move(curr_lengths);
        this->tau_inv = std::move(curr_tau_inv);
        this->max_length_ = new_max_length;
        this->intervals_ = this->lengths.size();
    }

    void apply_splitting(IntVectorType& curr_lengths, IntVectorType& curr_tau_inv, ulint& new_max_length) {
        if (this->split_params_ == NO_SPLITTING) {
            new_max_length = this->max_length_;
            return;
        }
        
        SplitResult<IntVectorType> split_result;
        
        if (this->split_params_.length_capping) {
            split_by_length_capping(curr_lengths, curr_tau_inv, this->domain_, *this->split_params_.length_capping, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_tau_inv = std::move(split_result.tau_inv);
            new_max_length = split_result.max_length;
        }
        if (this->split_params_.balancing) {
            split_by_balancing(curr_lengths, curr_tau_inv, this->domain_, *this->split_params_.balancing, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_tau_inv = std::move(split_result.tau_inv);
            new_max_length = split_result.max_length;
        }
    }
};

#endif // end include guard _PERMUTATION_HPP