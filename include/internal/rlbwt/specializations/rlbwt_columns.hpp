#ifndef _RLBWT_COLUMNS_HPP
#define _RLBWT_COLUMNS_HPP

#include "internal/common.hpp"
#include "internal/move/move_columns.hpp"

// Supports π = LF or π = FL, using new fields specific to RLBWT
// Need to extend the existing MoveCols enum to work properly with MoveColsTraits
// Uses RELATIVE positions using length and offset
enum class RLBWTCols {
    LENGTH, // RLBWT Length
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

// Supports π = LF or π = FL, using new fields specific to RLBWT
// We need to extend the existing MoveColsIdx enum to work properly with MoveColsTraits
// Supports ABSOLUTE positions using start, pointer, and offset
enum class RLBWTColsIdx {
    START, // i is the start of the interval
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

template <>
struct MoveColsTraits<RLBWTCols> {
    static constexpr bool RELATIVE = true;
    static constexpr RLBWTCols PRIMARY = RLBWTCols::LENGTH;

    static constexpr RLBWTCols LENGTH = RLBWTCols::LENGTH;
    static constexpr RLBWTCols POINTER = RLBWTCols::POINTER;
    static constexpr RLBWTCols OFFSET = RLBWTCols::OFFSET;
    static constexpr RLBWTCols CHARACTER = RLBWTCols::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(RLBWTCols::COUNT);

    using Position = MovePosition<RELATIVE>::type;
};

template <>
struct MoveColsTraits<RLBWTColsIdx> {
    static constexpr bool RELATIVE = false;
    static constexpr RLBWTColsIdx PRIMARY = RLBWTColsIdx::START;

    static constexpr RLBWTColsIdx START = RLBWTColsIdx::START;
    static constexpr RLBWTColsIdx POINTER = RLBWTColsIdx::POINTER;
    static constexpr RLBWTColsIdx OFFSET = RLBWTColsIdx::OFFSET;
    static constexpr RLBWTColsIdx CHARACTER = RLBWTColsIdx::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(RLBWTColsIdx::COUNT);
    
    using Position = MovePosition<RELATIVE>::type;
};

/* They just extend the existing switch mechanism */
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