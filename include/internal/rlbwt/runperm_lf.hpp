#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "internal/ds/alphabet.hpp"

template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermLFImpl : public RunPermRLBWT<RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, TableType>,
                         RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, TableType> {
    using Base = RunPermRLBWT<RunPermLFImpl, RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType, TableType>;
    using BaseColumns = typename Base::BaseColumns;
public:
    using Base::Base;
    using Base::operator=;

    void find_permutation_and_alphabet(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        AlphabetType& alphabet,
        ulint& num_chars,
        PackedVector<BaseColumns>& base_structure
    ) {
        auto [char_count, head_ranks, bwt_length] = get_char_counts(rlbwt_heads, rlbwt_run_lengths);
        num_chars = bwt_length;
        std::vector<ulint> interval_permutation = get_rlbwt_head_permutations(rlbwt_heads, char_count, head_ranks);
        alphabet = AlphabetType(char_count);
        auto mapped_rlbwt_heads = alphabet.map_sequence(rlbwt_heads);
        base_structure = Base::MoveStructureBase::find_structure(mapped_rlbwt_heads, rlbwt_run_lengths, interval_permutation, num_chars, alphabet.size());
    }

    void LF() {
        Base::next();
    }

    void LF(ulint steps) {
        Base::next(steps);
    }

private:
    // === Constructor Helpers ===
    static std::tuple<std::vector<size_t>, std::vector<size_t>, ulint> get_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
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

    std::vector<ulint> get_rlbwt_head_permutations(const std::vector<uchar> &rlbwt_heads, const std::vector<size_t> &char_count, const std::vector<size_t> &head_ranks) {
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
};

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MoveLFImpl : public MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, AlphabetType, TableType>, 
                      StoreAbsolutePositions, AlphabetType, TableType> {
    using Base = MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, AlphabetType, TableType>,
                 StoreAbsolutePositions, AlphabetType, TableType>;
public:
    using Base::Base;
    using Base::operator=;

    void LF() {
        Base::next();
    }

    void LF(ulint steps) {
        Base::next(steps);
    }
};

#endif /* end of include guard: _RLBWT_LF_HPP */