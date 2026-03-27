#ifndef _PHI_HELPERS_HPP
#define _PHI_HELPERS_HPP

#include "orbit/common.hpp"
#include "orbit/internal/ds/packed_vector_aligned.hpp"
#include "orbit/internal/rlbwt/lf_permutation.hpp"

namespace orbit::rlbwt {

using int_vec = int_vector_aligned;

template<typename lf_t>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_images(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, lf_t& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    int_vec phi_lengths(lf.runs(), bit_width(lf.domain() - 1));

    size_t UNUSED_INTERVAL = max_val(bit_width(lf.intervals()));
    size_t UNUSED_SA = max_val(bit_width(lf.domain()));
    int_vec move_run_to_phi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its Phi interval (only set those corresponding to RLBWT runs)
    int_vec run_tail_sa_samples(lf.intervals(), bit_width(lf.domain())); // The SA samples at the tail of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // Phi intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_phi_interval = lf.runs() - 1;
    // Step through entire BWT to recover Phi structure and SA samples at tails
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT runhead
        if (offset == 0) {
            if (interval == 0 || lf.get_character(interval - 1) != lf.get_character(interval)) {
                phi_lengths[curr_phi_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(phi_lengths[curr_phi_interval]));
                move_run_to_phi.set(interval, curr_phi_interval);
                last_sample = sa;
                --curr_phi_interval;
            }
            else {
                move_run_to_phi.set(interval, UNUSED_INTERVAL);
            }
        }
        // If at BWT run tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                run_tail_sa_samples.set(interval, sa);
            }
            else {
                run_tail_sa_samples.set(interval, UNUSED_SA);
            }
        }
        --sa;
        pos = lf.LF(pos);
    }

    int_vec phi_images(lf.runs(), bit_width(lf.domain() - 1));
    // Step through BWT tail samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_tail_sa_samples.size(); ++i) {
        if (run_tail_sa_samples.get(i) == UNUSED_SA) continue;
        phi_images[move_run_to_phi.get((i == lf.intervals() - 1) ? 0 : i + 1)] = run_tail_sa_samples.get(i);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_lengths, phi_images};
}

template<typename alphabet_t=nucleotide>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_images(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    lf_move_impl_default<alphabet_t> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi_images(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}   

template<typename lf_t>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_img_rank_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, lf_t& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    int_vec phi_lengths(lf.runs(), bit_width(lf.domain() - 1));

    size_t UNUSED_INTERVAL = max_val(bit_width(lf.intervals()));
    size_t UNUSED_VISIT_RANK = max_val(bit_width(lf.runs()));
    int_vec move_run_to_phi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its Phi interval (only set those corresponding to RLBWT runs)
    int_vec run_tail_visit_rank(lf.intervals(), bit_width(lf.runs())); // The rank of the visit to the tail of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // Phi intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_phi_interval = lf.runs() - 1;
    size_t curr_visit_rank = lf.runs() - 1;
    // Step through entire BWT to recover Phi structure and SA samples at tails
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT runhead
        if (offset == 0) {
            if (interval == 0 || lf.get_character(interval - 1) != lf.get_character(interval)) {
                phi_lengths[curr_phi_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(phi_lengths[curr_phi_interval]));
                move_run_to_phi.set(interval, curr_phi_interval);
                last_sample = sa;
                --curr_phi_interval;
            }
            else {
                move_run_to_phi.set(interval, UNUSED_INTERVAL);
            }
        }
        // If at BWT run tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                run_tail_visit_rank.set(interval, curr_visit_rank);
                --curr_visit_rank;
            }
            else {
                run_tail_visit_rank.set(interval, UNUSED_VISIT_RANK);
            }
        }
        --sa;
        pos = lf.LF(pos);
    }

    int_vec phi_img_rank_inv(lf.runs(), bit_width(lf.runs() - 1));
    // Step through BWT tail samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_tail_visit_rank.size(); ++i) {
        if (run_tail_visit_rank.get(i) == UNUSED_VISIT_RANK) continue;
        phi_img_rank_inv[run_tail_visit_rank.get(i)] = move_run_to_phi.get((i == lf.intervals() - 1) ? 0 : i + 1);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_lengths, phi_img_rank_inv};
}

template<typename alphabet_t=nucleotide>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_img_rank_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    lf_move_impl_default<alphabet_t> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi_img_rank_inv(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

template<typename lf_t>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_inv_images(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, lf_t& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    int_vec phi_inv_lengths(lf.runs(), bit_width(lf.domain() - 1));
    
    size_t UNUSED_INTERVAL = max_val(bit_width(lf.intervals()));
    size_t UNUSED_SA = max_val(bit_width(lf.domain()));
    int_vec move_run_to_phi_inv(lf.intervals(), bit_width(lf.intervals())); // Map move run to its phi_inv interval (only set those corresponding to RLBWT runs)
    int_vec run_head_sa_samples(lf.intervals(), bit_width(lf.domain())); // The SA samples at the head of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // phi_inv intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_phi_inv_interval = lf.runs() - 1;
    // Step through entire BWT to recover phi_inv structure and SA samples at heads
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                phi_inv_lengths[curr_phi_inv_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(phi_inv_lengths[curr_phi_inv_interval]));
                move_run_to_phi_inv.set(interval, curr_phi_inv_interval);
                last_sample = sa;
                --curr_phi_inv_interval;
            }
            else {
                move_run_to_phi_inv.set(interval, UNUSED_INTERVAL);
            }
        }
        // If at BWT run head
        if (offset == 0) {
            if (interval == 0 || lf.get_character(interval - 1) != lf.get_character(interval)) {
                run_head_sa_samples.set(interval, sa);
            }
            else {
                run_head_sa_samples.set(interval, UNUSED_SA);
            }
        }
        --sa;
        pos = lf.LF(pos);
    }

    int_vec phi_inv_images(lf.runs(), bit_width(lf.domain() - 1));
    // Step through BWT head samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_head_sa_samples.size(); ++i) {
        if (run_head_sa_samples.get(i) == UNUSED_SA) continue;
        phi_inv_images[move_run_to_phi_inv.get((i == 0) ? lf.intervals() - 1 : i - 1)] = run_head_sa_samples.get(i);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_inv_lengths, phi_inv_images};
}

template<typename alphabet_t=nucleotide>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_inv_images(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    lf_move_impl_default<alphabet_t> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi_inv_images(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

template<typename lf_t>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_inv_img_rank_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, lf_t& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    int_vec phi_inv_lengths(lf.runs(), bit_width(lf.domain() - 1));
    
    size_t UNUSED_INTERVAL = max_val(bit_width(lf.intervals()));
    size_t UNUSED_VISIT_RANK = max_val(bit_width(lf.runs()));
    int_vec move_run_to_phi_inv(lf.intervals(), bit_width(lf.intervals())); // Map move run to its phi_inv interval (only set those corresponding to RLBWT runs)
    int_vec run_head_visit_rank(lf.intervals(), bit_width(lf.runs())); // The rank of the visit to the head of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // phi_inv intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_phi_inv_interval = lf.runs() - 1;
    size_t curr_visit_rank = lf.runs() - 1;
    // Step through entire BWT to recover phi_inv structure and SA samples at heads
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                phi_inv_lengths[curr_phi_inv_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(phi_inv_lengths[curr_phi_inv_interval]));
                move_run_to_phi_inv.set(interval, curr_phi_inv_interval);
                last_sample = sa;
                --curr_phi_inv_interval;
            }
            else {
                move_run_to_phi_inv.set(interval, UNUSED_INTERVAL);
            }
        }
        // If at BWT run head
        if (offset == 0) {
            if (interval == 0 || lf.get_character(interval - 1) != lf.get_character(interval)) {
                run_head_visit_rank.set(interval, curr_visit_rank);
                --curr_visit_rank;
            }
            else {
                run_head_visit_rank.set(interval, UNUSED_VISIT_RANK);
            }
        }
        --sa;
        pos = lf.LF(pos);
    }

    int_vec phi_inv_img_rank_inv(lf.runs(), bit_width(lf.runs() - 1));
    // Step through BWT head samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_head_visit_rank.size(); ++i) {
        if (run_head_visit_rank.get(i) == UNUSED_VISIT_RANK) continue;
        phi_inv_img_rank_inv[run_head_visit_rank.get(i)] = move_run_to_phi_inv.get((i == 0) ? lf.intervals() - 1 : i - 1);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_inv_lengths, phi_inv_img_rank_inv};
}

template<typename alphabet_t=nucleotide>
inline std::tuple<int_vec, int_vec> rlbwt_to_phi_inv_img_rank_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    lf_move_impl_default<alphabet_t> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi_inv_img_rank_inv(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

} // namespace orbit::rlbwt

#endif /* end of include guard: _PHI_HELPERS_HPP */