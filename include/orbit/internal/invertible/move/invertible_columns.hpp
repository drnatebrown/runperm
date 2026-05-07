#ifndef ORBIT_INVERTIBLE_COLUMNS_HPP
#define ORBIT_INVERTIBLE_COLUMNS_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/move_columns.hpp"

namespace orbit {

enum class invertible_columns {
    LENGTH, // Length of move interval
    POINTER_FWD, // Pointer to forward move interval
    POINTER_INV, // Pointer to inverse move interval
    // Do not need offset for invertible move structure
    FWD_INTERVAL, // Interval was originally a forward move interval
    INV_INTERVAL, // Interval was originally an inverse move interval
    COUNT // Helper to get the number of columns
};

enum class invertible_columns_idx {
    START, // i is the start of the interval
    POINTER_FWD, // Pointer to forward move interval
    POINTER_INV, // Pointer to inverse move interval
    // Do not need offset for invertible move structure
    FWD_INTERVAL, // Interval was originally a forward move interval
    INV_INTERVAL, // Interval was originally an inverse move interval
    COUNT // Helper to get the number of columns
};

template<>
struct move_cols_traits<invertible::invertible_columns> {
    static constexpr bool RELATIVE = true;
    static constexpr invertible::invertible_columns PRIMARY = invertible::invertible_columns::LENGTH;

    static constexpr invertible::invertible_columns LENGTH = invertible::invertible_columns::LENGTH;
    static constexpr invertible::invertible_columns POINTER_FWD = invertible::invertible_columns::POINTER_FWD;
    static constexpr invertible::invertible_columns POINTER_INV = invertible::invertible_columns::POINTER_INV;
    static constexpr invertible::invertible_columns FWD_INTERVAL = invertible::invertible_columns::FWD_INTERVAL;
    static constexpr invertible::invertible_columns INV_INTERVAL = invertible::invertible_columns::INV_INTERVAL;
    static constexpr size_t NUM_COLS = static_cast<size_t>(invertible::invertible_columns::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct move_cols_traits<invertible::invertible_columns_idx> {
    static constexpr bool RELATIVE = false;
    static constexpr invertible::invertible_columns_idx PRIMARY = invertible::invertible_columns_idx::START;

    static constexpr invertible::invertible_columns_idx START = invertible::invertible_columns_idx::START;
    static constexpr invertible::invertible_columns_idx POINTER_FWD = invertible::invertible_columns_idx::POINTER_FWD;
    static constexpr invertible::invertible_columns_idx POINTER_INV = invertible::invertible_columns_idx::POINTER_INV;
    static constexpr invertible::invertible_columns_idx FWD_INTERVAL = invertible::invertible_columns_idx::FWD_INTERVAL;
    static constexpr invertible::invertible_columns_idx INV_INTERVAL = invertible::invertible_columns_idx::INV_INTERVAL;
    static constexpr size_t NUM_COLS = static_cast<size_t>(invertible::invertible_columns_idx::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct column_switcher<invertible::invertible_columns> {
    using relative = invertible::invertible_columns;
    using absolute = invertible::invertible_columns_idx;
};

template<>
struct column_switcher<invertible::invertible_columns_idx> {
    using relative = invertible::invertible_columns;
    using absolute = invertible::invertible_columns_idx;
};

} // namespace orbit
#endif /* end of include guard: ORBIT_INVERTIBLE_COLUMNS_HPP */