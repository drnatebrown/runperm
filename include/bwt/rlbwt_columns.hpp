#ifndef _RLBWT_COLUMNS_HPP
#define _RLBWT_COLUMNS_HPP

#include "common.hpp"
#include "move/move_columns.hpp"

// Supports π = LF or π = FL, using new fields specific to RLBWT
// Need to extend the existing MoveCols enum to work properly with MoveColsTraits
// Uses RELATIVE positions using length and offset
enum class RLBWTCols {
    LENGTH, // RLBWT Length
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    NUM_COLS // Helper to get the number of columns
};

// Supports π = LF or π = FL, using new fields specific to RLBWT
// We need to extend the existing MoveColsIdx enum to work properly with MoveColsTraits
// Supports ABSOLUTE positions using start, pointer, and offset
enum class RLBWTColsIdx {
    START, // i is the start of the interval
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    NUM_COLS // Helper to get the number of columns
};

/* They just extend the existing MoveCols enum */
template<>
struct MoveColsTraits<RLBWTCols> : public ExtendTraits<RLBWTCols, MoveCols> {};
template<>
struct MoveColsTraits<RLBWTColsIdx> : public ExtendTraits<RLBWTColsIdx, MoveColsIdx> {};

/* Similarly, they just extend the existing switch mechanism */
template<>
struct ColumnSwitcher<RLBWTCols> {
    using Relative = RLBWTCols;
    using Absolute = RLBWTColsIdx;
};

template<>
struct ColumnSwitcher<RLBWTColsIdx> {
    using Relative = RLBWTCols;
    using Absolute = RLBWTColsIdx;
};

#endif /* end of include guard: _RLBWT_COLUMNS_HPP */