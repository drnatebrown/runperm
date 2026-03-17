#ifndef _RLBWT_HELPERS_HPP
#define _RLBWT_HELPERS_HPP

#include "common.hpp"

namespace lf {

    inline std::tuple<std::vector<size_t>, std::vector<size_t>, ulint> get_LF_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());
    
        std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
        std::vector<size_t> head_ranks(rlbwt_heads.size(), 0);
        size_t bwt_length = 0;
        for (size_t i = 0; i < rlbwt_heads.size(); i++)
        {
            uchar c = rlbwt_heads[i];
            ulint length = rlbwt_run_lengths[i];
            head_ranks[i] = char_count[c];
            char_count[c] += length;
            bwt_length+=length;
        }
        return {char_count, head_ranks, bwt_length};
    }
    
    inline std::tuple<std::vector<size_t>, size_t, ulint> get_LF_head_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());
    
        std::vector<size_t> head_counts(MAX_ALPHABET_SIZE, 0);
        ulint max_length = 0;
        size_t bwt_length = 0;
        for (size_t i = 0; i < rlbwt_heads.size(); i++)
        {
            uchar c = rlbwt_heads[i];
            ulint length = rlbwt_run_lengths[i];
            ++head_counts[c];
            bwt_length+=length;
            max_length = std::max(max_length, length);
        }
        return {head_counts, bwt_length, max_length};
    }

    inline std::vector<ulint> get_LF_head_permutations(const std::vector<uchar> &rlbwt_heads, const std::vector<size_t> &char_count, const std::vector<size_t> &head_ranks) {
        std::vector<ulint> interval_permutation(head_ranks.size());
        
        std::vector<size_t> C_array(char_count.size(), 0);
        size_t seen = 0;
        for (size_t i = 0; i < char_count.size(); i++) {
            C_array[i] = seen;
            seen += char_count[i];
        }

        for (size_t i = 0; i < rlbwt_heads.size(); i++) {
            uchar c = rlbwt_heads[i];
            interval_permutation[i] = C_array[c] + head_ranks[i];
        }
        return interval_permutation;
    }
    
    inline std::vector<ulint> get_LF_tau_inv(const std::vector<uchar> &rlbwt_heads, const std::vector<size_t> &head_counts) {
        // TODO use IntVector for tau_inv
        std::vector<ulint> tau_inv(rlbwt_heads.size());
        std::vector<size_t> seen_heads(head_counts.size(), 0);
    
        std::vector<size_t> C_head_array(head_counts.size(), 0);
        size_t seen = 0;
        for (size_t i = 0; i < head_counts.size(); i++) {
            C_head_array[i] = seen;
            seen += head_counts[i];
        }
    
        for (size_t i = 0; i < rlbwt_heads.size(); i++) {
            uchar c = rlbwt_heads[i];
            tau_inv[C_head_array[c] + seen_heads[c]] = i;
            ++seen_heads[c];
        }
    
        return tau_inv;
    }
    
} // end namespace lf

namespace fl {
    std::tuple<std::vector<size_t>, std::vector<std::vector<std::pair<size_t, size_t>>>, ulint> get_FL_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
        // Holds the lengths of runs in F, and the origins of the runs in BWT (i.e., FL(i) for i at run heads)
        std::vector<std::vector<std::pair<size_t, size_t>>> F_lens_and_origins(MAX_ALPHABET_SIZE);
        size_t bwt_length = 0;
        for (size_t i = 0; i < rlbwt_heads.size(); i++)
        {
            uchar c = rlbwt_heads[i];
            ulint length = rlbwt_run_lengths[i];
            F_lens_and_origins[c].push_back({length, bwt_length});
            char_count[c] += length;
            bwt_length+=length;
        }
        return {char_count, F_lens_and_origins, bwt_length};
    }

    std::tuple<std::vector<size_t>, std::vector<std::vector<std::pair<size_t, size_t>>>, size_t, ulint> get_FL_head_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        std::vector<size_t> head_counts(MAX_ALPHABET_SIZE, 0);
        // Holds the lengths of runs in F, and the origin run in BWT (i.e., run containing the FL(i) for i at run heads)
        std::vector<std::vector<std::pair<size_t, size_t>>> F_lens_and_origin_run(MAX_ALPHABET_SIZE);
        size_t bwt_length = 0;
        ulint max_length = 0;
        for (size_t i = 0; i < rlbwt_heads.size(); i++)
        {
            uchar c = rlbwt_heads[i];
            ulint length = rlbwt_run_lengths[i];
            F_lens_and_origin_run[c].push_back({length, i});
            ++head_counts[c];
            bwt_length+=length;
            max_length = std::max(max_length, length);
        }
        return {head_counts, F_lens_and_origin_run, bwt_length, max_length};
    }

    std::tuple<std::vector<uchar>, std::vector<ulint>, std::vector<ulint>> get_FL_runs_and_interval_permutation(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins) {
        std::vector<uchar> F_heads(runs);
        std::vector<ulint> F_lens(runs);
        std::vector<ulint> interval_permutation(runs);

        size_t curr_run = 0;
        for (size_t c = 0; c < F_lens_and_origins.size(); ++c) {
            for (size_t j = 0; j < F_lens_and_origins[c].size(); ++j) {
                F_heads[curr_run] = c;
                F_lens[curr_run] = F_lens_and_origins[c][j].first;
                interval_permutation[curr_run] = F_lens_and_origins[c][j].second;
                curr_run++;
            }
        }
        return {F_heads, F_lens, interval_permutation};
    }

    std::tuple<std::vector<uchar>, std::vector<ulint>, std::vector<ulint>> get_FL_runs_and_tau_inv(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins) {
        std::vector<uchar> F_heads(runs);
        std::vector<ulint> F_lens(runs);
        std::vector<ulint> F_tau_inv(runs);

        size_t curr_run = 0;
        for (size_t c = 0; c < F_lens_and_origins.size(); ++c) {
            for (size_t j = 0; j < F_lens_and_origins[c].size(); ++j) {
                F_heads[curr_run] = c;
                F_lens[curr_run] = F_lens_and_origins[c][j].first;
                F_tau_inv[F_lens_and_origins[c][j].second] = curr_run;
                curr_run++;
            }
        }
        return {F_heads, F_lens, F_tau_inv};
    }
    
} // end namespace fl

#endif