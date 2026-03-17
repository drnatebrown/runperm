#ifndef _RUNPERM_INVPHI_HPP
#define _RUNPERM_INVPHI_HPP

#include "common.hpp"
#include "internal/runperm/runperm.hpp"
#include "internal/rlbwt/runperm_lf.hpp"

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


    MoveInvPhi(const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        const SplitParams& split_params = SplitParams())
    : Base([&] {
        size_t domain = 0;
        ulint max_length = 0;
        auto [invphi_lengths, invphi_tau_inv] =
            invphi::rlbwt_to_invphi_tau_inv<>(rlbwt_heads, rlbwt_run_lengths,
                                            &domain, &max_length);
        return PermutationImpl<>::from_lengths_and_tau_inv(
            invphi_lengths, invphi_tau_inv, domain, max_length, split_params);
    }()) {}

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