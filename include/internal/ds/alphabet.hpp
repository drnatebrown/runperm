#ifndef _ALPHABET_MAP_HPP
#define _ALPHABET_MAP_HPP

#include "internal/common.hpp"
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

    std::vector<uchar> map_sequence(const std::vector<uchar>& sequence) {
        std::vector<uchar> mapped_sequence(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
            mapped_sequence[i] = map_char(sequence[i]);
        }
        return mapped_sequence;
    }   

    std::vector<uchar> unmap_sequence(const std::vector<uchar>& sequence) {
        std::vector<uchar> unmapped_sequence(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
            unmapped_sequence[i] = unmap_char(sequence[i]);
        }
        return unmapped_sequence;
    }

    size_t serialize(std::ostream& out) {
        size_t written_bytes = 0;

        size_t size = alphabet_map.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        char* data = reinterpret_cast<char*>(alphabet_map.data());
        size = alphabet_map.size() * sizeof(uchar);
        out.write(data, size);
        written_bytes += size;

        size = reverse_alphabet_map.size();
        out.write((char *)&size, sizeof(size));
        written_bytes += sizeof(size);

        data = reinterpret_cast<char*>(reverse_alphabet_map.data());
        size = reverse_alphabet_map.size() * sizeof(uchar);
        out.write(data, size);
        written_bytes += size;

        return written_bytes;
    }

    void load(std::istream& in) {
        size_t size;
        in.read((char *)&size, sizeof(size));
        alphabet_map = std::vector<uchar>(size);
        char* data = reinterpret_cast<char*>(alphabet_map.data());
        size_t bytes = alphabet_map.size() * sizeof(uchar);
        in.read(data, bytes);

        in.read((char *)&size, sizeof(size));
        reverse_alphabet_map = std::vector<uchar>(size);
        data = reinterpret_cast<char*>(reverse_alphabet_map.data());
        bytes = reverse_alphabet_map.size() * sizeof(uchar);
        in.read(data, bytes);
    }

private:
    std::vector<uchar> alphabet_map;
    std::vector<uchar> reverse_alphabet_map;
};

class Nucleotide {
    public:
    static constexpr uchar SIGMA = 7;

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

    static std::vector<uchar> map_sequence(const std::vector<uchar>& sequence) {
        std::vector<uchar> mapped_sequence(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
            mapped_sequence[i] = map_char(sequence[i]);
        }
        return mapped_sequence;
    }   

    static std::vector<uchar> unmap_sequence(const std::vector<uchar>& sequence) {
        std::vector<uchar> unmapped_sequence(sequence.size());
        for (size_t i = 0; i < sequence.size(); i++) {
            unmapped_sequence[i] = unmap_char(sequence[i]);
        }
        return unmapped_sequence;
    }

    // Dummy methods, nothing to serialize or load
    size_t serialize(std::ostream& out) { return 0; }
    void load(std::istream& in) { }

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
        TER,SEP,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,NUL,  
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
    TER,SEP,'A','C','G','T','N'
    };
};

#endif