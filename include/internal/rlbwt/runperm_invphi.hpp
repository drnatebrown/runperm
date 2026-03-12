#ifndef _RUNPERM_INVPHI_HPP
#define _RUNPERM_INVPHI_HPP

#include "internal/common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/rlbwt/runperm_lf.hpp"

namespace phi {

using IntVec = IntVectorAligned;

template<typename LFType>
std::tuple<IntVec, IntVec> rlbwt_to_invphi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, LFType& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    IntVec invphi_lengths(lf.runs(), bit_width(lf.domain() - 1));
    
    size_t UNUSED_INTERVAL = MAX_VAL(bit_width(lf.intervals()));
    size_t UNUSED_SA = MAX_VAL(bit_width(lf.domain()));
    IntVec move_run_to_invphi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its InvPhi interval (only set those corresponding to RLBWT runs)
    IntVec run_head_sa_samples(lf.intervals(), bit_width(lf.domain())); // The SA samples at the head of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // InvPhi intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_invphi_interval = lf.runs() - 1;
    // Step through entire BWT to recover InvPhi structure and SA samples at heads
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                invphi_lengths[curr_invphi_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(invphi_lengths[curr_invphi_interval]));
                move_run_to_invphi.set(interval, curr_invphi_interval);
                last_sample = sa;
                --curr_invphi_interval;
            }
            else {
                move_run_to_invphi.set(interval, UNUSED_INTERVAL);
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

    IntVec invphi_interval_permutations(lf.runs(), bit_width(lf.domain() - 1));
    // Step through BWT head samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_head_sa_samples.size(); ++i) {
        if (run_head_sa_samples.get(i) == UNUSED_SA) continue;
        invphi_interval_permutations[move_run_to_invphi.get((i == 0) ? lf.intervals() - 1 : i - 1)] = run_head_sa_samples.get(i);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {invphi_lengths, invphi_interval_permutations};
}

template<typename AlphabetType=Nucleotide>
inline std::tuple<IntVec, IntVec> rlbwt_to_invphi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    MoveLFImplDefault<AlphabetType> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_invphi(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

template<typename LFType>
std::tuple<IntVec, IntVec> rlbwt_to_invphi_tau_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, LFType& lf, size_t* domain = nullptr, ulint* max_length = nullptr) {
    IntVec invphi_lengths(lf.runs(), bit_width(lf.domain() - 1));
    
    size_t UNUSED_INTERVAL = MAX_VAL(bit_width(lf.intervals()));
    size_t UNUSED_VISIT_RANK = MAX_VAL(bit_width(lf.runs()));
    IntVec move_run_to_invphi(lf.intervals(), bit_width(lf.intervals())); // Map move run to its InvPhi interval (only set those corresponding to RLBWT runs)
    IntVec run_head_visit_rank(lf.intervals(), bit_width(lf.runs())); // The rank of the visit to the head of each move run (only set those corresponding to RLBWT runs)

    ulint max_length_seen = 0;
    auto pos = lf.first();
    size_t last_sample = lf.domain();
    size_t sa = lf.domain() - 1;
    // InvPhi intervals correspond to the original (unsplit) permutation runs, not move runs.
    size_t curr_invphi_interval = lf.runs() - 1;
    size_t curr_visit_rank = lf.runs() - 1;
    // Step through entire BWT to recover InvPhi structure and SA samples at heads
    for (size_t i = 0; i < lf.domain(); ++i) {
        size_t interval = pos.interval;
        size_t offset = pos.offset;
        // If at BWT tail
        if (offset == lf.get_length(interval) - 1) {
            if (interval == lf.intervals() - 1 || lf.get_character(interval + 1) != lf.get_character(interval)) {
                invphi_lengths[curr_invphi_interval] = last_sample - sa;
                max_length_seen = std::max(max_length_seen, static_cast<ulint>(invphi_lengths[curr_invphi_interval]));
                move_run_to_invphi.set(interval, curr_invphi_interval);
                last_sample = sa;
                --curr_invphi_interval;
            }
            else {
                move_run_to_invphi.set(interval, UNUSED_INTERVAL);
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

    IntVec invphi_tau_inv(lf.runs(), bit_width(lf.runs() - 1));
    // Step through BWT head samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_head_visit_rank.size(); ++i) {
        if (run_head_visit_rank.get(i) == UNUSED_VISIT_RANK) continue;
        invphi_tau_inv[run_head_visit_rank.get(i)] = move_run_to_invphi.get((i == 0) ? lf.intervals() - 1 : i - 1);
    }

    if (domain != nullptr) {
        *domain = lf.domain();
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }

    return {invphi_lengths, invphi_tau_inv};
}

template<typename AlphabetType=Nucleotide>
inline std::tuple<IntVec, IntVec> rlbwt_to_invphi_tau_inv(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    // Need a move structure with LF to find SA samples
    MoveLFImplDefault<AlphabetType> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_invphi_tau_inv(bwt_heads, bwt_run_lengths, move_lf, domain, max_length);
}

} // end namespace invphi

// Always use absolute positions for InvPhi
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         template<typename> class TableType = MoveVector>
class RunPermInvPhi : public RunPermImpl<RunColsType, IntegratedMoveStructure, true, ExponentialSearch, MoveCols, MoveStructure, TableType> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, true, ExponentialSearch, MoveCols, MoveStructure, TableType>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position InvPhi(Position pos) {
        return Base::next(pos);
    }

    Position InvPhi(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

    ulint SA(Position pos) {
        return pos.idx;
    }
};

// Always use absolute positions for InvPhi
class MoveInvPhi : public MovePermImpl<true, DEFAULT_EXPONENTIAL_SEARCH, MoveCols, MoveStructure, MoveVector> {
    using Base = MovePermImpl<true, DEFAULT_EXPONENTIAL_SEARCH, MoveCols, MoveStructure, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position InvPhi(Position pos) {
        return Base::next(pos);
    }

    Position InvPhi(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

    ulint SA(Position pos) {
        return pos.idx;
    }
};


#endif /* end of include guard: _RUNPERM_INVPHI_HPP */