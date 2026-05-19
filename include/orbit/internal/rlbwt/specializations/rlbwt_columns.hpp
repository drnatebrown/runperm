#ifndef _RLBWT_COLUMNS_HPP
#define _RLBWT_COLUMNS_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/move_columns.hpp"

namespace orbit::rlbwt {

// Supports π = LF or π = FL, using new fields specific to RLBWT
// Need to extend the existing move_columns enum to work properly with move_cols_traits
// Uses RELATIVE positions using length and offset
enum class rlbwt_columns {
    LENGTH, // RLBWT Length
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

// Supports π = LF or π = FL, using new fields specific to RLBWT
// We need to extend the existing move_columns_idx enum to work properly with move_cols_traits
// Supports ABSOLUTE positions using start, pointer, and offset
enum class rlbwt_columns_idx {
    START, // i is the start of the interval
    POINTER, // Corresponds to character that is BWT run head, interval is the run containing its π
    OFFSET, // Offset of π within its run
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

} // namespace orbit::rlbwt

namespace orbit {

template <>
struct move_cols_traits<rlbwt::rlbwt_columns> {
    static constexpr bool RELATIVE = true;
    static constexpr bool INVERTIBLE = false;
    static constexpr rlbwt::rlbwt_columns PRIMARY = rlbwt::rlbwt_columns::LENGTH;

    static constexpr rlbwt::rlbwt_columns LENGTH = rlbwt::rlbwt_columns::LENGTH;
    static constexpr rlbwt::rlbwt_columns POINTER = rlbwt::rlbwt_columns::POINTER;
    static constexpr rlbwt::rlbwt_columns OFFSET = rlbwt::rlbwt_columns::OFFSET;
    static constexpr rlbwt::rlbwt_columns CHARACTER = rlbwt::rlbwt_columns::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(rlbwt::rlbwt_columns::COUNT);

    using position = move_position<RELATIVE>::type;
};

template <>
struct move_cols_traits<rlbwt::rlbwt_columns_idx> {
    static constexpr bool RELATIVE = false;
    static constexpr bool INVERTIBLE = false;
    static constexpr rlbwt::rlbwt_columns_idx PRIMARY = rlbwt::rlbwt_columns_idx::START;

    static constexpr rlbwt::rlbwt_columns_idx START = rlbwt::rlbwt_columns_idx::START;
    static constexpr rlbwt::rlbwt_columns_idx POINTER = rlbwt::rlbwt_columns_idx::POINTER;
    static constexpr rlbwt::rlbwt_columns_idx OFFSET = rlbwt::rlbwt_columns_idx::OFFSET;
    static constexpr rlbwt::rlbwt_columns_idx CHARACTER = rlbwt::rlbwt_columns_idx::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(rlbwt::rlbwt_columns_idx::COUNT);
    
    using position = move_position<RELATIVE>::type;
};

/* They just extend the existing switch mechanism */
template<>
struct column_switcher<rlbwt::rlbwt_columns> {
    using relative = rlbwt::rlbwt_columns;
    using absolute = rlbwt::rlbwt_columns_idx;
};

template<>
struct column_switcher<rlbwt::rlbwt_columns_idx> {
    using relative = rlbwt::rlbwt_columns;
    using absolute = rlbwt::rlbwt_columns_idx;
};

} // namespace orbit

// ================================ RLBWT INVERTIBLE MOVE STRUCTURES ================================

namespace orbit::rlbwt {

enum class rlbwt_invertible_columns {
    LENGTH, // Length of move interval
    POINTER_FWD, // Pointer to forward move interval
    POINTER_INV, // Pointer to inverse move interval
    // Do not need offset for invertible move structure
    FWD_INTERVAL, // Interval was originally a forward move interval
    INV_INTERVAL, // Interval was originally an inverse move interval
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

enum class rlbwt_invertible_columns_idx {
    START, // i is the start of the interval
    POINTER_FWD, // Pointer to forward move interval
    POINTER_INV, // Pointer to inverse move interval
    // Do not need offset for invertible move structure
    FWD_INTERVAL, // Interval was originally a forward move interval
    INV_INTERVAL, // Interval was originally an inverse move interval
    CHARACTER, // RLBWT Character
    COUNT // Helper to get the number of columns
};

} // namespace orbit::rlbwt

namespace orbit {

template<>
struct move_cols_traits<rlbwt::rlbwt_invertible_columns> {
    static constexpr bool RELATIVE = true;
    static constexpr bool INVERTIBLE = true;
    static constexpr rlbwt::rlbwt_invertible_columns PRIMARY = rlbwt::rlbwt_invertible_columns::LENGTH;

    static constexpr rlbwt::rlbwt_invertible_columns LENGTH = rlbwt::rlbwt_invertible_columns::LENGTH;
    static constexpr rlbwt::rlbwt_invertible_columns POINTER_FWD = rlbwt::rlbwt_invertible_columns::POINTER_FWD;
    static constexpr rlbwt::rlbwt_invertible_columns POINTER_INV = rlbwt::rlbwt_invertible_columns::POINTER_INV;
    static constexpr rlbwt::rlbwt_invertible_columns FWD_INTERVAL = rlbwt::rlbwt_invertible_columns::FWD_INTERVAL;
    static constexpr rlbwt::rlbwt_invertible_columns INV_INTERVAL = rlbwt::rlbwt_invertible_columns::INV_INTERVAL;
    static constexpr rlbwt::rlbwt_invertible_columns CHARACTER = rlbwt::rlbwt_invertible_columns::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(rlbwt::rlbwt_invertible_columns::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct move_cols_traits<rlbwt::rlbwt_invertible_columns_idx> {
    static constexpr bool RELATIVE = false;
    static constexpr bool INVERTIBLE = true;
    static constexpr rlbwt::rlbwt_invertible_columns_idx PRIMARY = rlbwt::rlbwt_invertible_columns_idx::START;

    static constexpr rlbwt::rlbwt_invertible_columns_idx START = rlbwt::rlbwt_invertible_columns_idx::START;
    static constexpr rlbwt::rlbwt_invertible_columns_idx POINTER_FWD = rlbwt::rlbwt_invertible_columns_idx::POINTER_FWD;
    static constexpr rlbwt::rlbwt_invertible_columns_idx POINTER_INV = rlbwt::rlbwt_invertible_columns_idx::POINTER_INV;
    static constexpr rlbwt::rlbwt_invertible_columns_idx FWD_INTERVAL = rlbwt::rlbwt_invertible_columns_idx::FWD_INTERVAL;
    static constexpr rlbwt::rlbwt_invertible_columns_idx INV_INTERVAL = rlbwt::rlbwt_invertible_columns_idx::INV_INTERVAL;
    static constexpr rlbwt::rlbwt_invertible_columns_idx CHARACTER = rlbwt::rlbwt_invertible_columns_idx::CHARACTER;
    static constexpr size_t NUM_COLS = static_cast<size_t>(rlbwt::rlbwt_invertible_columns_idx::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct column_switcher<rlbwt::rlbwt_invertible_columns> {
    using relative = rlbwt::rlbwt_invertible_columns;
    using absolute = rlbwt::rlbwt_invertible_columns_idx;
};

template<>
struct column_switcher<rlbwt::rlbwt_invertible_columns_idx> {
    using relative = rlbwt::rlbwt_invertible_columns;
    using absolute = rlbwt::rlbwt_invertible_columns_idx;
};

} // namespace orbit

#endif /* end of include guard: _RLBWT_COLUMNS_HPP */