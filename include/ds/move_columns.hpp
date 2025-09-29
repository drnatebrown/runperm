#ifndef _MOVE_COLUMNS_HPP
#define _MOVE_COLUMNS_HPP

#include "common.hpp"

/* Stores permutation π(1..n) in O(r)-space, where r is the number of positions such that
either i == 0 or π(i-1) != π(i) - 1. An interval of length l is any range where π(i) is stored, and
otherwise π(j), for i <=j < (i - l), is computed as π(j) = π(i) + (j - i). Here, we use interval/offset notation. */
enum class MoveCols { 
    LENGTH, // Length of move interval
    POINTER, // Where i is start of current interval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    NUM_COLS // Helper to get the number of columns
};

/* Above, but uses absolute indices instead of relative intervals */
enum class MoveColsIdx { 
    START, // i is the start of the interval
    POINTER, // Where i is start of current inverval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    NUM_COLS // Helper to get the number of columns
};

template <typename ColumnsType>
struct MoveColsTraits;

template <>
struct MoveColsTraits<MoveCols> {
    static constexpr bool IS_LENGTH = true;
    static constexpr bool IS_START = !IS_LENGTH;
    static constexpr auto PRIMARY = MoveCols::LENGTH;
    static constexpr size_t NUM_COLS = static_cast<size_t>(MoveCols::NUM_COLS);

    struct Position {
        ulint interval = 0;
        ulint offset = 0;
    };
};

template <>
struct MoveColsTraits<MoveColsIdx> {
    static constexpr bool IS_LENGTH = false;
    static constexpr bool IS_START = true;
    static constexpr auto PRIMARY = MoveColsIdx::START;
    static constexpr size_t NUM_COLS = static_cast<size_t>(MoveColsIdx::NUM_COLS);
    
    struct Position {
        ulint interval = 0;
        ulint offset = 0;
        ulint idx = 0;
    };
};

#endif /* end of include guard: _MOVE_COLUMNS_HPP */