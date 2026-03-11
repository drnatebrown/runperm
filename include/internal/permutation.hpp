#ifndef _PERMUTATION_HPP
#define _PERMUTATION_HPP

#include "internal/common.hpp"
#include "internal/ds/packed_vector.hpp"
#include "internal/ds/packed_vector_aligned.hpp"
#include "internal/move/move_splitting.hpp"

#include <optional>
#include <numeric>
#include <algorithm>
#include <cassert>

namespace permutation_helpers {

inline std::tuple<std::vector<ulint>, ulint> starts_to_lengths(const std::vector<ulint>& starts, const ulint domain) {
    std::vector<ulint> lengths(starts.size());
    ulint max_length = 0;
    for (size_t i = 0; i < lengths.size() - 1; ++i) {
        lengths[i] = starts[i + 1] - starts[i];
        max_length = std::max(max_length, lengths[i]);
    }
    lengths[lengths.size() - 1] = domain - starts[starts.size() - 1];
    max_length = std::max(max_length, lengths[lengths.size() - 1]);
    return {lengths, max_length};
}

inline std::tuple<ulint, ulint> sum_and_max(const std::vector<ulint>& data) {
    ulint sum = 0;
    ulint max = 0;
    for (const auto& element : data) {
        sum += element;
        max = std::max(max, element);
    }
    return {sum, max};
}

inline std::vector<ulint> get_tau_inv(const std::vector<ulint>& interval_output_starts) {
    std::vector<ulint> tau_inv(interval_output_starts.size());
    std::iota(tau_inv.begin(), tau_inv.end(), 0);
    std::sort(tau_inv.begin(), tau_inv.end(), [&](ulint a, ulint b) { return interval_output_starts[a] < interval_output_starts[b]; });
    return tau_inv;
}

} // end namespace permutation_helpers

template<class IntVectorType = IntVectorAligned>
class PermutationImpl {
public:

    PermutationImpl() = default;

    static PermutationImpl<IntVectorType> from_permutation(const std::vector<ulint>& permutation, const SplitParams& split_params = SplitParams()) {
        ulint max_length = 0;
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation, &max_length);
        assert(lengths.size() == interval_permutation.size());

        return from_lengths_and_interval_permutation(lengths, interval_permutation, permutation.size(), max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_lengths_and_tau_inv(const std::vector<ulint>& lengths, const std::vector<ulint>& tau_inv, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau_inv.size());
        auto [domain, max_length] = permutation_helpers::sum_and_max(lengths);
        
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }
    
    static PermutationImpl<IntVectorType> from_lengths_and_tau_inv(const std::vector<ulint>& lengths, const std::vector<ulint>& tau_inv, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau_inv.size());
        PermutationImpl<IntVectorType> permutation;
        permutation.set_initial_values(domain, lengths.size(), max_length);
        permutation.init_tau_inv(lengths, tau_inv, split_params);
        return permutation;
    }

    static PermutationImpl<IntVectorType> from_lengths_and_tau(const std::vector<ulint>& lengths, const std::vector<ulint>& tau, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau.size());

        auto [domain, max_length] = permutation_helpers::sum_and_max(lengths);
        return from_lengths_and_tau(lengths, tau, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_lengths_and_tau(const std::vector<ulint>& lengths, const std::vector<ulint>& tau, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == tau.size());

        PermutationImpl<IntVectorType> permutation;
        permutation.set_initial_values(domain, lengths.size(), max_length);
        permutation.init_tau(lengths, tau, split_params);
        return permutation;
    }

    static PermutationImpl<IntVectorType> from_lengths_and_interval_permutation(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        auto [domain, max_length] = permutation_helpers::sum_and_max(lengths);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_lengths_and_interval_permutation(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        std::vector<ulint> tau_inv = permutation_helpers::get_tau_inv(interval_permutation);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_tau_inv(const std::vector<ulint>& starts, const std::vector<ulint>& tau_inv, const ulint domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau_inv.size());
        auto [lengths, max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_tau_inv(const std::vector<ulint>& starts, const std::vector<ulint>& tau_inv, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau_inv.size());
        auto [lengths, calculated_max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_tau_inv(lengths, tau_inv, domain, calculated_max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_tau(const std::vector<ulint>& starts, const std::vector<ulint>& tau, const ulint domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau.size());
        auto [lengths, max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        return from_lengths_and_tau(lengths, tau, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_tau(const std::vector<ulint>& starts, const std::vector<ulint>& tau, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == tau.size());
        auto [lengths, calculated_max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_tau(lengths, tau, domain, calculated_max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_interval_permutation(const std::vector<ulint>& starts, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == interval_permutation.size());
        auto [lengths, max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, max_length, split_params);
    }

    static PermutationImpl<IntVectorType> from_starts_and_interval_permutation(const std::vector<ulint>& starts, const std::vector<ulint>& interval_permutation, const ulint domain, const ulint max_length, const SplitParams& split_params = SplitParams()) {
        assert(starts.size() == interval_permutation.size());
        auto [lengths, calculated_max_length] = permutation_helpers::starts_to_lengths(starts, domain);
        assert(calculated_max_length == max_length);
        return from_lengths_and_interval_permutation(lengths, interval_permutation, domain, calculated_max_length, split_params);
    }

    template<typename T>
    std::vector<T> split_run_data_with_copy(const std::vector<ulint>& original_lengths, const std::vector<T>& run_data) {
        assert(original_lengths.size() == run_data.size());
        assert(original_lengths.size() == original_runs);

        std::vector<T> new_run_data(split_runs);
        size_t new_interval_idx = 0;
        for (size_t i = 0; i < original_lengths.size(); ++i) {
            ulint remaining = original_lengths[i];
            while (remaining > 0) {
                ulint chunk = this->lengths[new_interval_idx];
                new_run_data[new_interval_idx] = run_data[i];
                remaining -= chunk;
                ++new_interval_idx;
            }
        }
        assert(new_interval_idx == split_runs);

        return new_run_data;
    }

    ulint get_length(size_t i) {
        return lengths[i];
    }
    ulint get_tau_inv(size_t i) {
        return tau_inv[i];
    }
    size_t get_domain() {
        return domain;
    }
    size_t get_original_runs() {
        return original_runs;
    }
    size_t get_split_runs() {
        return split_runs;
    }
    ulint get_max_length() {
        return max_length;
    }

private:
    IntVectorType lengths;
    IntVectorType tau_inv;
    size_t domain;
    size_t original_runs;
    size_t split_runs;
    ulint max_length;

    void set_initial_values(size_t domain, size_t original_runs, ulint max_length) {
        this->domain = domain;
        this->original_runs = original_runs;
        this->max_length = max_length;
        this->split_runs = original_runs;
    }
    
    void init_tau(const std::vector<ulint>& lengths, const std::vector<ulint>& tau, const SplitParams& split_params) {
        uchar length_bits = bit_width(max_length);
        uchar tau_inv_bits = bit_width(original_runs - 1);

        IntVectorType curr_lengths(original_runs, length_bits);
        IntVectorType curr_tau_inv(original_runs, tau_inv_bits);
        for (size_t i = 0; i < original_runs; ++i) {
            curr_lengths[i] = lengths[i];
            curr_tau_inv[tau[i]] = i;
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_tau_inv, new_max_length, split_params);

        this->lengths = std::move(curr_lengths);
        this->tau_inv = std::move(curr_tau_inv);
        this->max_length = new_max_length;
        this->split_runs = this->lengths.size();
    }

    void init_tau_inv(const std::vector<ulint>& lengths, const std::vector<ulint>& tau_inv, const SplitParams& split_params) {
        uchar length_bits = bit_width(max_length);
        uchar tau_inv_bits = bit_width(original_runs - 1);

        IntVectorType curr_lengths(original_runs, length_bits);
        IntVectorType curr_tau_inv(original_runs, tau_inv_bits);
        for (size_t i = 0; i < original_runs; ++i) {
            curr_lengths[i] = lengths[i];
            curr_tau_inv[i] = tau_inv[i];
        }

        ulint new_max_length = 0;
        apply_splitting(curr_lengths, curr_tau_inv, new_max_length, split_params);

        this->lengths = std::move(curr_lengths);
        this->tau_inv = std::move(curr_tau_inv);
        this->max_length = new_max_length;
        this->split_runs = this->lengths.size();
    }

    void apply_splitting(IntVectorType& curr_lengths, IntVectorType& curr_tau_inv, ulint& new_max_length, const SplitParams& split_params) {
        if (!split_params.length_capping && !split_params.balancing) {
            new_max_length = max_length;
            return;
        }
        
        SplitResult<IntVectorType> split_result;
        
        if (split_params.length_capping) {
            split_by_length_capping(curr_lengths, curr_tau_inv, this->domain, *split_params.length_capping, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_tau_inv = std::move(split_result.tau_inv);
            new_max_length = split_result.max_length;
        }
        if (split_params.balancing) {
            split_by_balancing(curr_lengths, curr_tau_inv, domain, *split_params.balancing, split_result);
            curr_lengths = std::move(split_result.lengths);
            curr_tau_inv = std::move(split_result.tau_inv);
            new_max_length = split_result.max_length;
        }
    }
};

#endif // end include guard _PERMUTATION_HPP