#ifndef _RUNPERM_PHI_HPP
#define _RUNPERM_PHI_HPP

#include "internal/common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/rlbwt/runperm_lf.hpp"

namespace phi {

using IntVec = IntVectorAligned;

template<typename LFType>
std::tuple<IntVec, IntVec> rlbwt_to_phi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, LFType& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    IntVec phi_lengths(lf.runs(), bit_width(lf.domain() - 1));

    size_t UNUSED_INTERVAL = MAX_VAL(bit_width(lf.intervals()));
    size_t UNUSED_SA = MAX_VAL(bit_width(lf.domain()));
    IntVec move_run_to_phi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its Phi interval (only set those corresponding to RLBWT runs)
    IntVec run_tail_sa_samples(lf.intervals(), bit_width(lf.domain())); // The SA samples at the tail of each move run (only set those corresponding to RLBWT runs)

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

    IntVec phi_interval_permutations(lf.runs(), bit_width(lf.domain() - 1));
    // Step through BWT tail samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_tail_sa_samples.size(); ++i) {
        if (run_tail_sa_samples.get(i) == UNUSED_SA) continue;
        phi_interval_permutations[move_run_to_phi.get((i == lf.intervals() - 1) ? 0 : i + 1)] = run_tail_sa_samples.get(i);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_lengths, phi_interval_permutations};
}

template<typename AlphabetType=Nucleotide>
inline std::tuple<IntVec, IntVec> rlbwt_to_phi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    MoveLFImplDefault<AlphabetType> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

template<typename LFType>
std::tuple<IntVec, IntVec> rlbwt_to_phi_tau_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, LFType& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    IntVec phi_lengths(lf.runs(), bit_width(lf.domain() - 1));

    size_t UNUSED_INTERVAL = MAX_VAL(bit_width(lf.intervals()));
    size_t UNUSED_VISIT_RANK = MAX_VAL(bit_width(lf.runs()));
    IntVec move_run_to_phi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its Phi interval (only set those corresponding to RLBWT runs)
    IntVec run_tail_visit_rank(lf.intervals(), bit_width(lf.runs())); // The rank of the visit to the tail of each move run (only set those corresponding to RLBWT runs)

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

    IntVec phi_tau_inv(lf.runs(), bit_width(lf.runs() - 1));
    // Step through BWT tail samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_tail_visit_rank.size(); ++i) {
        if (run_tail_visit_rank.get(i) == UNUSED_VISIT_RANK) continue;
        phi_tau_inv[run_tail_visit_rank.get(i)] = move_run_to_phi.get((i == lf.intervals() - 1) ? 0 : i + 1);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {phi_lengths, phi_tau_inv};
}

template<typename AlphabetType=Nucleotide>
inline std::tuple<IntVec, IntVec> rlbwt_to_phi_tau_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    MoveLFImplDefault<AlphabetType> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi_tau_inv(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

} // end namespace phi

// Always use absolute positions for Phi
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class TableType = MoveVector>
class RunPermPhi : public RunPermImpl<RunColsType, IntegratedMoveStructure, true, ExponentialSearch, MoveCols, MoveStructure, TableType> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, true, ExponentialSearch, MoveCols, MoveStructure, TableType>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position Phi(Position pos) {
        return Base::next(pos);
    }

    Position Phi(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

    ulint SA(Position pos) {
        return pos.idx;
    }
};

// Always use absolute positions for Phi
class MovePhi : public MovePermImpl<true, DEFAULT_EXPONENTIAL_SEARCH, MoveCols, MoveStructure, MoveVector> {
    using Base = MovePermImpl<true, DEFAULT_EXPONENTIAL_SEARCH, MoveCols, MoveStructure, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    MovePhi(const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const SplitParams& split_params = SplitParams())
    : Base([&] {
        size_t domain = 0;
        ulint max_length = 0;
        auto [phi_lengths, phi_tau_inv] =
            phi::rlbwt_to_phi_tau_inv<>(rlbwt_heads, rlbwt_run_lengths,
                                        &domain, &max_length);
        return Permutation::from_lengths_and_tau_inv(
            phi_lengths, phi_tau_inv, domain, max_length, split_params);
      }()) {}

    Position Phi(Position pos) {
        return Base::next(pos);
    }

    Position Phi(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

    ulint SA(Position pos) {
        return pos.idx;
    }
};


#endif /* end of include guard: _RUNPERM_PHI_HPP */