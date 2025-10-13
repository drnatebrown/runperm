#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "internal/ds/alphabet.hpp"

template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename BaseColumnsType = RLBWTCols,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermFL : public RunPermRLBWT<RunPermFL<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, BaseColumnsType, AlphabetType, TableType>,
                         RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, BaseColumnsType, AlphabetType, TableType> {
    using Base = RunPermRLBWT<RunPermFL, RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, BaseColumnsType, AlphabetType, TableType>;
public:

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
        std::vector<uchar> alphabet_characters = Base::apply_alphabet(F_heads);
        base_structure = Base::MoveStructureBase::find_structure(alphabet_characters, F_lens, interval_permutation, num_chars, alphabet.size());
    }

    void FL() {
        Base::move();
    }

    void FL(ulint steps) {
        Base::move(steps);
    }

private:
    // === Constructor Helpers ===
    static std::tuple<std::vector<size_t>, std::vector<size_t>, ulint> get_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
        assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

        std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
        // Holds the lengths of runs in F, and the origins of the runs in BWT (i.e., FL(i) for i at run heads)
        std::vector<std::vector<std::pair<size_t, size_t>>> F_lens_and_origins(MAX_ALPHABET_SIZE);
        size_t bwt_length = 0;
        for (size_t i = 0; i < rlbwt_heads.size(); i++)
        {
            uchar c = rlbwt_heads[i];
            ulint length = rlbwt_run_lengths[i];
            if (c <= TERMINATOR) c = TERMINATOR;
            else if (c > TERMINATOR && c <= SEPARATOR) c = SEPARATOR;

            F_lens_and_origins[c].push_back({length, bwt_length});
            char_count[c] += length;
            bwt_length+=length;
        }
        return {char_count, F_lens_and_origins, bwt_length};
    }

    std::vector<ulint> get_F_runs(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins) {
        std::vector<uchar> F_heads(runs);
        std::vector<ulint> F_lens(runs);
        std::vector<ulint> interval_permutation(runs));

        size_t curr_run = 0;
        for (size_t i = 0; i < F_lens_and_origins.size(); i++) {
            for (size_t j = 0; j < F_lens_and_origins[i].size(); j++) {
                F_heads[curr_run] = j;
                F_lens[curr_run] = F_lens_and_origins[i][j].first;
                interval_permutation[curr_run] = F_lens_and_origins[i][j].second;
                curr_run++;
            }
        }
        return {F_heads, F_lens, interval_permutation};
    }
};

template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLIntegrated = RunPermFL<RunColsType, true, DEFAULT_STORE_ABSOLUTE_POSITIONS, RLBWTCols, AlphabetType, MoveVector>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLIntegratedAbsolute = RunPermFL<RunColsType, true, true, RLBWTCols, AlphabetType, MoveVector>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLSeperated = RunPermFL<RunColsType, false, DEFAULT_STORE_ABSOLUTE_POSITIONS, RLBWTCols, AlphabetType, MoveVector>;
template<typename RunColsType, typename AlphabetType = Nucleotide>
using RunPermFLSeperatedAbsolute = RunPermFL<RunColsType, false, true, RLBWTCols, AlphabetType, MoveVector>;

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename BaseColumns = RLBWTCols, 
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MoveFL : public MovePermRLBWT<RunPermFL, StoreAbsolutePositions, BaseColumns, AlphabetType, TableType> {
public:
    using Base = MovePermRLBWT<RunPermFL, StoreAbsolutePositions, BaseColumns, AlphabetType, TableType>;

    void FL() {
        Base::FL();
    }

    void FL(ulint steps) {
        Base::FL(steps);
    }
};

#endif /* end of include guard: _RLBWT_LF_HPP */