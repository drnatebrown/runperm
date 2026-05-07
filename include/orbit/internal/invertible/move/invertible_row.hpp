#ifndef _MOVE_ROW_HPP
#define _MOVE_ROW_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/invertible/invertible_columns.hpp"
#include "orbit/internal/move/move_row.hpp"
#include <stdexcept>
#include <array>
#include <cassert>

namespace orbit {

template <typename columns_t, size_t PRIMARY_BITS, size_t POINTER_FWD_BITS, size_t POINTER_INV_BITS>
struct invertible_row_bits {};

template <typename columns_t, size_t P, size_t PTR_FWD, size_t PTR_INV>
struct move_row_traits<invertible::invertible_row_bits<columns_t, P, PTR_FWD, PTR_INV>> {
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_FWD_BITS = PTR_FWD;
    static constexpr size_t POINTER_INV_BITS = PTR_INV;

    static constexpr size_t FWD_INTERVAL_BITS = 1;
    static constexpr size_t INV_INTERVAL_BITS = 1;
};

template <typename columns_t = move_columns>
struct move_row {
    // Sets NumCols, Columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)
    using row_traits = move_row_traits<columns>;

    ulint primary : row_traits::PRIMARY_BITS;
    ulint pointer_fwd : row_traits::POINTER_FWD_BITS;
    ulint pointer_inv : row_traits::POINTER_INV_BITS;
    ulint fwd_interval : row_traits::FWD_INTERVAL_BITS;
    ulint inv_interval : row_traits::INV_INTERVAL_BITS;

    move_row() = default;
    move_row(const std::array<ulint, num_cols>& values) {
        set(values);
    }

    template <columns col>
    void set(ulint val) {
        if constexpr (col == cols_traits::PRIMARY) primary = val;
        else if constexpr (col == cols_traits::POINTER_FWD) pointer_fwd = val;
        else if constexpr (col == cols_traits::POINTER_INV) pointer_inv = val;
        else if constexpr (col == cols_traits::FWD_INTERVAL) fwd_interval = val;
        else if constexpr (col == cols_traits::INV_INTERVAL) inv_interval = val;
    }
    template <size_t... indices>
    void set(const std::array<ulint, num_cols>& values, std::index_sequence<indices...>) {
        (set<static_cast<columns>(indices)>(values[indices]), ...);
    }
    void set(const std::array<ulint, num_cols>& values) {
        set(values, std::make_index_sequence<num_cols>{});
    }

    template <columns col>
    ulint get() const {
        if constexpr (col == cols_traits::PRIMARY) return primary;
        else if constexpr (col == cols_traits::POINTER_FWD) return pointer_fwd;
        else if constexpr (col == cols_traits::POINTER_INV) return pointer_inv;
        else if constexpr (col == cols_traits::FWD_INTERVAL) return fwd_interval;
        else if constexpr (col == cols_traits::INV_INTERVAL) return inv_interval;
        else throw std::invalid_argument("Invalid column");
    }
    template <size_t... indices>
    std::array<ulint, num_cols> get(std::index_sequence<indices...>) const {
        return {get<static_cast<columns>(indices)>()...};
    }
    std::array<ulint, num_cols> get() const {
        return get(std::make_index_sequence<num_cols>{});
    }

    static const std::array<uchar, num_cols>& get_widths() {
        static const std::array<uchar, num_cols> widths{
            static_cast<uchar>(row_traits::PRIMARY_BITS),
            static_cast<uchar>(row_traits::POINTER_FWD_BITS),
            static_cast<uchar>(row_traits::POINTER_INV_BITS),
            static_cast<uchar>(row_traits::FWD_INTERVAL_BITS),
            static_cast<uchar>(row_traits::INV_INTERVAL_BITS),
        };
        return widths;
    }

    static void assert_widths(const std::array<uchar, num_cols>& widths) {
        assert(widths[static_cast<size_t>(cols_traits::PRIMARY)] <= row_traits::PRIMARY_BITS);
        assert(widths[static_cast<size_t>(cols_traits::POINTER_FWD)] <= row_traits::POINTER_FWD_BITS);
        assert(widths[static_cast<size_t>(cols_traits::POINTER_INV)] <= row_traits::POINTER_INV_BITS);
        assert(widths[static_cast<size_t>(cols_traits::FWD_INTERVAL)] <= row_traits::FWD_INTERVAL_BITS);
        assert(widths[static_cast<size_t>(cols_traits::INV_INTERVAL)] <= row_traits::INV_INTERVAL_BITS);
    }

} __attribute__((packed));

// Column Sizes for move_columns (+ 2 bits for fwd/inv interval)
using invertible_cols_4 = invertible_row_bits<invertible_columns, 8, 14, 8>;     // (4 bytes)
using invertible_cols_5 = invertible_row_bits<invertible_columns, 8, 22, 8>;     // (5 bytes)
using invertible_cols_6 = invertible_row_bits<invertible_columns, 12, 22, 12>; // (6 bytes)
using invertible_cols_7 = invertible_row_bits<invertible_columns, 12, 30, 12>; // (7 bytes)
using invertible_cols_8 = invertible_row_bits<invertible_columns, 16, 30, 16>; // (8 bytes)
using invertible_cols_9 = invertible_row_bits<invertible_columns, 16, 38, 16>; // (9 bytes)
using invertible_cols_10 = invertible_row_bits<invertible_columns, 20, 38, 20>; // (10 bytes)
using invertible_cols_11 = invertible_row_bits<invertible_columns, 20, 46, 20>; // (11 bytes)
using invertible_cols_12 = invertible_row_bits<invertible_columns, 24, 46, 24>; // (12 bytes)
using invertible_cols_default = invertible_cols_9;

// Column Sizes for move_columns_idx (+ 2 bits for fwd/inv interval)
using invertible_cols_idx_5 = invertible_row_bits<invertible_columns_idx, 16, 14, 8>; // (5 bytes)
using invertible_cols_idx_6 = invertible_row_bits<invertible_columns_idx, 20, 18, 8>; // (6 bytes)
using invertible_cols_idx_7 = invertible_row_bits<invertible_columns_idx, 26, 24, 8>; // (7 bytes)
using invertible_cols_idx_8 = invertible_row_bits<invertible_columns_idx, 28, 26, 12>; // (8 bytes)
using invertible_cols_idx_9 = invertible_row_bits<invertible_columns_idx, 32, 30, 12>; // (9 bytes)
using invertible_cols_idx_10 = invertible_row_bits<invertible_columns_idx, 34, 32, 12>; // (10 bytes)
using invertible_cols_idx_11 = invertible_row_bits<invertible_columns_idx, 36, 34, 16>; // (11 bytes)
using invertible_cols_idx_12 = invertible_row_bits<invertible_columns_idx, 40, 38, 16>; // (12 bytes)
using invertible_cols_idx_13 = invertible_row_bits<invertible_columns_idx, 44, 42, 16>; // (13 bytes)
using invertible_cols_idx_14 = invertible_row_bits<invertible_columns_idx, 46, 44, 20>; // (14 bytes)
using invertible_cols_idx_15 = invertible_row_bits<invertible_columns_idx, 50, 48, 20>; // (15 bytes)
using invertible_cols_idx_default = invertible_cols_idx_12;

template <>
struct move_row_traits<invertible_columns> : move_row_traits<invertible_cols_default> {};

template <>
struct move_row_traits<invertible_columns_idx> : move_row_traits<invertible_cols_idx_default> {};

} // namespace orbit

#endif /* end of include guard: _MOVE_ROW_HPP */