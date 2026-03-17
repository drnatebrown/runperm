#ifndef _RLBWT_PERMUTATION_HPP
#define _RLBWT_PERMUTATION_HPP

#include "common.hpp"
#include "internal/move/permutation_impl.hpp"
#include "internal/ds/alphabet.hpp"
#include "internal/rlbwt/rlbwt_helpers.hpp"

template<typename IntVectorType = IntVectorAligned, typename AlphabetType = Nucleotide>
class RLBWTPermutationImpl : public PermutationImpl<IntVectorType> {
    using Base = PermutationImpl<IntVectorType>;
public:
    using AlphabetTag = AlphabetType;

    RLBWTPermutationImpl() = default;

    static RLBWTPermutationImpl lf_permutation(const std::vector<uchar>& rlbwt_heads, const std::vector<ulint>& rlbwt_run_lengths, const SplitParams& split_params = SplitParams()) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        RLBWTPermutationImpl permutation;

        auto [head_counts, n, max_length] = lf::get_LF_head_counts(rlbwt_heads, rlbwt_run_lengths);
        permutation.set_initial_values(n, rlbwt_heads.size(), max_length, split_params);
        std::vector<ulint> tau_inv = lf::get_LF_tau_inv(rlbwt_heads, head_counts);

        permutation.init_tau_inv(rlbwt_run_lengths, tau_inv);

        permutation.alphabet_ = AlphabetType(head_counts);
        permutation.init_heads(rlbwt_heads, rlbwt_run_lengths);
        return permutation;
    }

    static RLBWTPermutationImpl fl_permutation(const std::vector<uchar>& rlbwt_heads, const std::vector<ulint>& rlbwt_run_lengths, const SplitParams& split_params = SplitParams()) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        RLBWTPermutationImpl permutation;

        auto [head_counts, F_lens_and_origin_run, n, max_length] = fl::get_FL_head_counts(rlbwt_heads, rlbwt_run_lengths);
        permutation.set_initial_values(n, rlbwt_heads.size(), max_length, split_params);

        auto [F_heads, F_lens, F_tau_inv] = fl::get_FL_runs_and_tau_inv(rlbwt_heads.size(), F_lens_and_origin_run);
        permutation.init_tau_inv(F_lens, F_tau_inv);

        permutation.alphabet_ = AlphabetType(head_counts);
        permutation.init_heads(F_heads, F_lens);
        return permutation;
    }

    const IntVectorType& get_heads() const {
        return heads_;
    }

    const AlphabetType& get_alphabet() const {
        return alphabet_;
    }

    const ulint sigma() const {
        return alphabet_.size();
    }

private:
    AlphabetType alphabet_;
    IntVectorType heads_;

    template<typename Collection1, typename Collection2>
    void init_heads(const Collection1& heads, const Collection2& original_lengths) {
        heads_ = IntVectorType(Base::intervals(), bit_width(alphabet_.size() - 1));
        size_t new_interval_idx = 0;
        for (size_t i = 0; i < original_lengths.size(); ++i) {
            ulint remaining = original_lengths[i];
            while (remaining > 0) {
                ulint chunk = Base::get_length(new_interval_idx);
                heads_.set(new_interval_idx, alphabet_.map_char(heads[i]));
                remaining -= chunk;
                ++new_interval_idx;
            }
        }
        assert(new_interval_idx == Base::intervals());
    }
};

#endif