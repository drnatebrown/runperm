#ifndef _RLBWT_PERMS_HPP
#define _RLBWT_PERMS_HPP

#include "common.hpp"
#include "runperm/runperm.hpp"
#include "rlbwt/alphabet.hpp"

enum class ColumnLF {
    CHARACTER,
    COUNT
};

template<typename AlphabetType = Nucleotide, typename RunPermType = RunPerm<ColumnLF>>
class MoveLF : public RunPermType {
public: 
    using RunData = typename RunPermType::RunData;

    MoveLF(const std::vector<uchar> &bwt_heads, const std::vector<ulint> &bwt_run_lengths, SplitParams split_params = SplitParams())
    {
        auto [char_count, head_ranks, num_chars] = get_char_counts(bwt_heads, bwt_run_lengths);
        std::vector<ulint> interval_permutation = get_bwt_head_permutations(bwt_heads, char_count, head_ranks);

        alphabet = AlphabetType(char_count);
        std::vector<RunData> alphabet_characters = apply_alphabet(bwt_heads);

        // Initialize base class
        if (split_params.max_allowed_length || split_params.balancing_factor) {
            static_cast<RunPermType&>(*this) = RunPermType(bwt_run_lengths, interval_permutation, num_chars, split_params, alphabet_characters);
        } else {
            static_cast<RunPermType&>(*this) = RunPermType(bwt_run_lengths, interval_permutation, num_chars, alphabet_characters);
        }
    }

    uchar get_character() {
        return alphabet.unmap_char(RunPermType::template get<ColumnLF::CHARACTER>());
    }

private:
    AlphabetType alphabet;
    
    // === Constructor Helpers ===
    static std::tuple<std::vector<size_t>, std::vector<size_t>, ulint> get_char_counts(const std::vector<uchar> &bwt_heads, const std::vector<ulint> &bwt_run_lengths) {
        assert(bwt_heads.size() == bwt_run_lengths.size());

        std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
        std::vector<size_t> head_ranks(bwt_heads.size(), 0);
        size_t num_chars = 0;
        for (size_t i = 0; i < bwt_heads.size(); i++)
        {
            uchar c = bwt_heads[i];
            ulint length = bwt_run_lengths[i];
            if (c <= TERMINATOR) c = TERMINATOR;
            else if (c > TERMINATOR && c <= SEPARATOR) c = SEPARATOR;
            head_ranks[i] = char_count[c];
            char_count[c] += length;
            num_chars+=length;
        }
        return {char_count, head_ranks, num_chars};
    }
    
    std::vector<RunData> apply_alphabet(const std::vector<uchar> &characters) {
        std::vector<RunData> transformed_characters(characters.size());
        for (size_t i = 0; i < characters.size(); i++) {
            transformed_characters[i][to_index(ColumnLF::CHARACTER)] = alphabet.map_char(characters[i]);
        }
        return transformed_characters;
    }

    std::vector<ulint> get_bwt_head_permutations(const std::vector<uchar> &bwt_heads, const std::vector<size_t> &char_count, const std::vector<size_t> &head_ranks) {
        std::vector<ulint> interval_permutation(head_ranks.size());
        
        
        std::vector<size_t> C_array(char_count.size(), 0);
        size_t seen = 0;
        for (size_t i = 0; i < char_count.size(); i++) {
            C_array[i] = seen;
            seen += char_count[i];
        }

        for (size_t i = 0; i < bwt_heads.size(); i++) {
            uchar c = bwt_heads[i];
            interval_permutation[i] = C_array[c] + head_ranks[i];
        }
        return interval_permutation;
    }
};

template<typename AlphabetType = Nucleotide, typename RunPermType = RunPerm<ColumnLF>>
class MoveFL : public RunPermType {
public: 
    using RunData = typename RunPermType::RunData;

    MoveFL(const std::vector<uchar> &bwt_heads, const std::vector<ulint> &bwt_run_lengths, SplitParams split_params = SplitParams())
    {
        
    }

private:
AlphabetType alphabet;
};
#endif /* end of include guard: _RLBWT_PERMS_HPP */