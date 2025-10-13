#ifndef _RUNPERM_PHI_HPP
#define _RUNPERM_PHI_HPP

#include "internal/common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/rlbwt/runperm_lf.hpp"

template<typename LFType>
std::tuple<std::vector<ulint>, std::vector<ulint>, size_t> rlbwt_to_phi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths, LFType& lf) {
    assert(bwt_heads.size() == bwt_run_lengths.size());

    // The lengths and interval permutations of the phi structure
    std::vector<ulint> phi_lengths(lf.permutation_runs());
    std::vector<ulint> phi_interval_permutations(lf.permutation_runs());

    // Helpers for computation
    IntVector move_run_to_phi(lf.move_runs(), bit_width(lf.move_runs())); // Map move run to its Phi interval (only set those corresponding to RLBWT runs)
    IntVector run_tail_sa_samples(lf.move_runs(), bit_width(lf.size())); // The SA samples at the tail of each move run (only set those corresponding to RLBWT runs)

    lf.first();
    size_t last_sample = lf.size();
    size_t sa = lf.size() - 1;
    size_t curr_phi_interval = lf.move_runs() - 1;
    // Step through entire BWT to recover Phi structure and SA samples at tails
    for (size_t i = 0; i < lf.size(); ++i) {
        size_t interval = lf.get_position().interval;
        size_t offset = lf.get_position().offset;
        // If at BWT runhead
        if (offset == 0 && (interval == 0 || lf.get_character(interval - 1) != lf.get_character(interval))) {
            phi_lengths[curr_phi_interval] = last_sample - sa;
            move_run_to_phi.set(interval, curr_phi_interval);
            last_sample = sa;
            --curr_phi_interval;
        }
        // If at BWT run tail
        if (offset == lf.get_length(interval) - 1 && (interval == lf.move_runs() - 1 || lf.get_character(interval + 1) != lf.get_character(interval))) {
            run_tail_sa_samples.set(interval, sa);
        }
        --sa;
        lf.LF();
    }

    // Step through BWT tail samples to fill in Phi interval permutations
    for (size_t i = 0; i < run_tail_sa_samples.size(); ++i) {
        phi_interval_permutations[move_run_to_phi.get((i == lf.move_runs() - 1) ? 0 : i + 1)] = run_tail_sa_samples.get(i);
    }

    return {phi_lengths, phi_interval_permutations, lf.size()};
}

std::tuple<std::vector<ulint>, std::vector<ulint>, size_t> rlbwt_to_phi(const std::vector<uchar>& bwt_heads, const std::vector<ulint>& bwt_run_lengths) {
    // Need a move structure with LF to find SA samples
    MoveLFImpl<> move_lf(bwt_heads, bwt_run_lengths);
    return rlbwt_to_phi(bwt_heads, bwt_run_lengths, move_lf);
}

// Always use absolute positions for Phi
template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         template<typename> class TableType = MoveVector>
class RunPermPhi : public RunPermImpl<RunColsType, IntegratedMoveStructure, true, MoveCols, MoveStructure, TableType> {
    using Base = RunPermImpl<RunColsType, IntegratedMoveStructure, true, MoveCols, MoveStructure, TableType>;
public:
    using Base::Base;
    using Base::operator=;

    void Phi() {
        Base::next();
    }

    void Phi(ulint steps) {
        Base::next(steps);
    }

    ulint SA() {
        return Base::get_position().idx;
    }
};

// Always use absolute positions for Phi
class MovePhi : public MovePermImpl<true, MoveCols, MoveStructure, MoveVector> {
    using Base = MovePermImpl<true, MoveCols, MoveStructure, MoveVector>;
public:
    using Base::Base;
    using Base::operator=;

    void Phi() {
        Base::next();
    }

    void Phi(ulint steps) {
        Base::next(steps);
    }

    ulint SA() {
        return Base::get_position().idx;
    }
};


#endif /* end of include guard: _RUNPERM_PHI_HPP */