#ifndef _COMMON_HPP
#define _COMMON_HPP
#include <stddef.h>
#include <utility>
#include <vector>

// LEAVE UNCHANGED
#define MAX_ALPHABET_SIZE 256
#define RW_BYTES 5

// CONFIGURABLE
#define TERMINATOR 0
#define SEPARATOR 1

#define START_BYTES 5
#define LENGTH_BYTES 2
#define POINTER_BYTES 4
#define OFFSET_BYTES 2
#define CHARACTER_BYTES 1
#define DEFAULT_BYTES 4

#define BYTES_TO_BITS(bytes) ((bytes) * 8)
#define BITS_TO_BYTES(bits) ((bits) / 8)
#define NUM_BITS(type) (BYTES_TO_BITS(sizeof(type))) 
#define CEIL_DIV(num, den) ((num + den - 1) / den)
#define POW2(bits) (1ULL << (bits))
#define MAX_VAL(bits) (POW2(bits) - 1)
#define MASK(bits) MAX_VAL(bits)

#define MOVE_STRUCTURE_EXTENSION ".move"

typedef unsigned char uchar;
typedef unsigned long int ulint;

constexpr uchar bit_width(ulint value) {
    return value == 0 ? 1 : 64 - __builtin_clzll(value);
}

// ENUM REPRESENTS COLUMNS, USE ENUM HELPERS TO ENFORCE STRUCTURE
template<class E>
constexpr size_t to_index(E e) noexcept { return static_cast<size_t>(e); }

#define MOVE_CLASS_TRAITS(ColumnsParam) \
    using Columns = ColumnsParam; \
    template<typename C> \
    using ColsTraitsFor = typename ResolveColsTraits<C>::type;  \
    using ColsTraits = typename ResolveColsTraits<Columns>::type;  \
    static constexpr size_t NumCols = ColsTraits::NUM_COLS; \
    template<typename E> static constexpr Columns to_cols(E e) { return static_cast<Columns>(e); }

/* PERMUTATION UTILITIES */
inline std::pair<std::vector<ulint>, std::vector<ulint>> get_permutation_intervals(const std::vector<ulint> &permutation) {
    std::vector<ulint> lengths;
    std::vector<ulint> interval_permutation;
    for (size_t i = 0; i < permutation.size(); ++i) {
        if (i == 0 || permutation[i] != permutation[i - 1] + 1) {
            lengths.push_back(1);
            interval_permutation.push_back(permutation[i]);
        } else {
            ++lengths.back();
        }
    }
    return {lengths, interval_permutation};
}

inline std::pair<std::vector<uchar>, std::vector<ulint>> bwt_to_rlbwt(const std::vector<uchar> &bwt_chars) {
    std::vector<uchar> rlbwt_chars;
    std::vector<ulint> rlbwt_run_lengths;
    for (size_t i = 0; i < bwt_chars.size(); ++i) {
        if (i == 0 || bwt_chars[i] != bwt_chars[i - 1]) {
            rlbwt_chars.push_back(bwt_chars[i]);
            rlbwt_run_lengths.push_back(1);
        } else {
            ++rlbwt_run_lengths.back();
        }
    }
    return {rlbwt_chars, rlbwt_run_lengths};
}

#endif /* end of include guard: _COMMON_HPP */