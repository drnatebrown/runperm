#ifndef _ALPHABET_MAP_HPP
#define _ALPHABET_MAP_HPP

#include "common.hpp"
#include <cctype>
#include <cassert>

static constexpr uchar UNMAPPED = MAX_VAL(BYTES_TO_BITS(sizeof(uchar)));

class Alphabet
{
public:

    Alphabet() = default;

    Alphabet(const std::vector<ulint> &char_counts) {
        assert(char_counts.size() == MAX_ALPHABET_SIZE);

        alphabet_map.resize(MAX_ALPHABET_SIZE, UNMAPPED);
        uchar map_value = 0;
        for (uchar c = 0; c < MAX_ALPHABET_SIZE; c++) {
            if (char_counts[c] > 0) {
                alphabet_map[c] = map_value;
                reverse_alphabet_map.push_back(c);
                map_value++;
            }
        }

        if (reverse_alphabet_map.size() == MAX_ALPHABET_SIZE) {
            throw std::runtime_error("Alphabet mapping is over entire max alphabet size, do not map!");
        }
    }

    uchar map_char(uchar c) {
        assert(alphabet_map[c] != UNMAPPED);
        return alphabet_map[c];
    }

    uchar unmap_char(uchar c) {
        assert(c < reverse_alphabet_map.size());
        return reverse_alphabet_map[c];
    }

    uchar size() {
        return reverse_alphabet_map.size();
    }

private:
    std::vector<uchar> alphabet_map;
    std::vector<uchar> reverse_alphabet_map;
};

class Nucleotide {
    public:
    static constexpr uchar SIGMA = 8;

    Nucleotide() = default;

    Nucleotide(const std::vector<ulint> &char_counts) {
        for (size_t i = 0; i < char_counts.size(); i++) {
            if (char_counts[i] > 0) {
                if(map_char(i) == NUL) {
                    throw std::runtime_error("Nucletotide mapping contains non-DNA character: " + std::to_string(i));
                }
            }
        }
    }

    static constexpr uchar map_char(uchar c) {
        assert(alphabet_map[c] != NUL);
        return alphabet_map[c];
    }

    static constexpr uchar unmap_char(uchar c) {
        assert(c < size());
        return reverse_alphabet_map[c];
    }

    static constexpr uchar size() {
        return SIGMA;
    }

    private:
    static constexpr uchar TER = 0;
    static constexpr uchar SEP = 1;
    static constexpr uchar A = 2;
    static constexpr uchar C = 3;
    static constexpr uchar G = 4;
    static constexpr uchar T = 5;
    static constexpr uchar N = 6;
    static constexpr uchar NUL = UNMAPPED;

    static constexpr uchar alphabet_map[MAX_ALPHABET_SIZE] = {
        NUL,TER,SEP,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL, 
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL, A ,NUL, C ,NUL,NUL,NUL, G ,NUL,NUL,NUL,NUL,NUL,NUL, N ,NUL, 
        NUL,NUL,NUL,NUL, T ,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
        NUL, A ,NUL, C ,NUL,NUL,NUL, G ,NUL,NUL,NUL,NUL,NUL,NUL, N ,NUL,  
        NUL,NUL,NUL,NUL, T ,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,
        NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL
    };

    static constexpr uchar reverse_alphabet_map[SIGMA] = {
    TER,SEP,'A','C','G','T','N',NUL
    };
};

#endif