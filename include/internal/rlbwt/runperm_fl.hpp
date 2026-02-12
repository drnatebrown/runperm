#ifndef _RLBWT_FL_HPP
#define _RLBWT_FL_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "internal/ds/alphabet.hpp"

template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermFLImpl : public RunPermRLBWT<RunPermFLImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
                         RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = RunPermRLBWT<RunPermFLImpl, RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
    using BaseColumns = typename Base::BaseColumns;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    void find_permutation_and_alphabet(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        AlphabetType& alphabet,
        ulint& num_chars,
        PackedVector<BaseColumns>& base_structure
    ) {
        auto [char_count, F_lens_and_origins, bwt_length] = get_char_counts(rlbwt_heads, rlbwt_run_lengths);
        num_chars = bwt_length;
        auto [F_heads, F_lens, interval_permutation] = get_F_runs(rlbwt_heads.size(), F_lens_and_origins);
        alphabet = AlphabetType(char_count);
        auto mapped_F_heads = alphabet.map_sequence(F_heads);
        base_structure = Base::MoveStructureBase::find_structure(mapped_F_heads, F_lens, interval_permutation, num_chars, alphabet.size());
    }

    Position FL(Position pos) {
        return Base::next(pos);
    }

    Position FL(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

private:
    // === Constructor Helpers ===
    static std::tuple<std::vector<size_t>, std::vector<std::vector<std::pair<size_t, size_t>>>, ulint> get_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
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

    std::tuple<std::vector<uchar>, std::vector<ulint>, std::vector<ulint>> get_F_runs(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins) {
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
};

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MoveFLImpl : public MovePermRLBWT<RunPermFLImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>, 
                StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = MovePermRLBWT<RunPermFLImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
    StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position FL(Position pos) {
        return Base::next(pos);
    }

    Position FL(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }
};

#endif /* end of include guard: _RLBWT_FL_HPP */