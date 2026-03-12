#ifndef _RLBWT_LF_HPP
#define _RLBWT_LF_HPP

#include "internal/common.hpp"
#include "internal/permutation.hpp"
#include "internal/rlbwt/specializations/runperm_rlbwt.hpp"
#include "internal/ds/alphabet.hpp"

// TODO use IntVector here

template<typename RunColsType,
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE,
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class RunPermLFImpl : public RunPermRLBWT<RunPermLFImpl<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
                         RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = RunPermRLBWT<RunPermLFImpl, RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
    using BaseColumns = typename Base::BaseColumns;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    // TODO use IntVector and container templates here
    void find_permutation_and_alphabet(
        const std::vector<uchar>& rlbwt_heads,
        const std::vector<ulint>& rlbwt_run_lengths,
        AlphabetType& alphabet,
        size_t& bwt_length,
        PackedVector<BaseColumns>& base_structure,
        const SplitParams& split_params
    ) {
        auto [head_counts, n, max_length] = get_head_counts(rlbwt_heads, rlbwt_run_lengths);
        std::vector<ulint> tau_inv = get_rlbwt_tau_inv(rlbwt_heads, head_counts);
        bwt_length = n;
        alphabet = AlphabetType(head_counts);

        std::vector<uchar> mapped_rlbwt_heads = alphabet.map_sequence(rlbwt_heads);
        Permutation permutation = Permutation::from_lengths_and_tau_inv(rlbwt_run_lengths, tau_inv, bwt_length, max_length, split_params);
        std::vector<uchar> split_rlbwt_heads = permutation.split_run_data_with_copy(rlbwt_run_lengths, mapped_rlbwt_heads);

        base_structure = Base::MoveStructureBase::find_structure(split_rlbwt_heads, permutation, alphabet.size());
    }

    Position LF(Position pos) {
        return Base::next(pos);
    }

    Position LF(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }

private:
    // === Constructor Helpers ===
    // static std::tuple<std::vector<size_t>, std::vector<size_t>, ulint> get_char_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
    //     assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

    //     std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
    //     std::vector<size_t> head_ranks(rlbwt_heads.size(), 0);
    //     size_t bwt_length = 0;
    //     for (size_t i = 0; i < rlbwt_heads.size(); i++)
    //     {
    //         uchar c = rlbwt_heads[i];
    //         ulint length = rlbwt_run_lengths[i];
    //         head_ranks[i] = char_count[c];
    //         char_count[c] += length;
    //         bwt_length+=length;
    //     }
    //     return {char_count, head_ranks, bwt_length};
    // }

    static std::tuple<std::vector<size_t>, size_t, ulint> get_head_counts(const std::vector<uchar> &rlbwt_heads, const std::vector<ulint> &rlbwt_run_lengths) {
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

    // std::vector<ulint> get_rlbwt_head_permutations(const std::vector<uchar> &rlbwt_heads, const std::vector<size_t> &char_count, const std::vector<size_t> &head_ranks) {
    //     std::vector<ulint> interval_permutation(head_ranks.size());
        
    //     std::vector<size_t> C_array(char_count.size(), 0);
    //     size_t seen = 0;
    //     for (size_t i = 0; i < char_count.size(); i++) {
    //         C_array[i] = seen;
    //         seen += char_count[i];
    //     }

    //     for (size_t i = 0; i < rlbwt_heads.size(); i++) {
    //         uchar c = rlbwt_heads[i];
    //         interval_permutation[i] = C_array[c] + head_ranks[i];
    //     }
    //     return interval_permutation;
    // }

    std::vector<ulint> get_rlbwt_tau_inv(const std::vector<uchar> &rlbwt_heads, const std::vector<size_t> &head_counts) {
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
};

template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS,
         bool ExponentialSearch = DEFAULT_EXPONENTIAL_SEARCH,
         typename AlphabetType = Nucleotide,
         template<typename> class TableType = MoveVector>
class MoveLFImpl : public MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>, 
                      StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType> {
    using Base = MovePermRLBWT<RunPermLFImpl<EmptyRunCols, false, StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>,
                 StoreAbsolutePositions, ExponentialSearch, AlphabetType, TableType>;
public:
    using Base::Base;
    using Base::operator=;
    using Position = typename Base::Position;

    Position LF(Position pos) {
        return Base::next(pos);
    }

    Position LF(Position pos, ulint steps) {
        return Base::next(pos, steps);
    }
};

template<typename Alphabet=Nucleotide>
using MoveLFImplDefault = MoveLFImpl<DEFAULT_STORE_ABSOLUTE_POSITIONS, DEFAULT_EXPONENTIAL_SEARCH, Alphabet, MoveVector>;

#endif /* end of include guard: _RLBWT_LF_HPP */