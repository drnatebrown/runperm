#ifndef _RUNPERM_PHI_HPP
#define _RUNPERM_PHI_HPP

#include "common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/rlbwt/phi_helpers.hpp"

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