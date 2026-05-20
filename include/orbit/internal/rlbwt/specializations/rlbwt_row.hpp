#ifndef _RLBWT_ROW_HPP
#define _RLBWT_ROW_HPP

#include "orbit/common.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "orbit/internal/move/move_row.hpp"
#include <stdexcept>
#include <array>

namespace orbit::rlbwt {

// rlbwt_row_bits takes bit values directly
template <typename columns_t, size_t CHARACTER_BITS, size_t PRIMARY_BITS, size_t POINTER_BITS, size_t OFFSET_BITS>
struct rlbwt_row_bits {
};

template <typename columns_t, size_t CHARACTER_BITS, size_t PRIMARY_BITS, size_t POINTER_BITS>
struct invertible_rlbwt_row_bits {};

} // namespace orbit::rlbwt

namespace orbit {

// Specialization for move_row_traits
template <typename columns_t, size_t C, size_t P, size_t PTR, size_t OFF>
struct move_row_traits<rlbwt::rlbwt_row_bits<columns_t, C, P, PTR, OFF>> {
    static constexpr size_t CHARACTER_BITS = C;
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_BITS = PTR;
    static constexpr size_t OFFSET_BITS = OFF;
};

} // namespace orbit

namespace orbit::rlbwt {

template <typename columns_t = rlbwt_columns>
struct rlbwt_row {
    // Sets num_cols, columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)
    using row_traits = move_row_traits<columns>;

    ulint primary : row_traits::PRIMARY_BITS;
    ulint pointer : row_traits::POINTER_BITS;
    ulint offset : row_traits::OFFSET_BITS;
    ulint character : row_traits::CHARACTER_BITS;

    rlbwt_row() = default;
    rlbwt_row(const std::array<ulint, num_cols>& values) {
        set(values);
    }

    template <columns col>
    void set(ulint val) {
        if constexpr (col == cols_traits::PRIMARY) primary = val;
        else if constexpr (col == cols_traits::POINTER) pointer = val;
        else if constexpr (col == cols_traits::OFFSET) offset = val;
        else if constexpr (col == cols_traits::CHARACTER) character = val;
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
        else if constexpr (col == cols_traits::POINTER) return pointer;
        else if constexpr (col == cols_traits::OFFSET)   return offset;
        else if constexpr (col == cols_traits::CHARACTER) return character;
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
            static_cast<uchar>(row_traits::CHARACTER_BITS),
            static_cast<uchar>(row_traits::PRIMARY_BITS),
            static_cast<uchar>(row_traits::POINTER_BITS),
            static_cast<uchar>(row_traits::OFFSET_BITS),
            static_cast<uchar>(row_traits::CHARACTER_BITS),
        };
        return widths;
    }

    static void assert_widths(const std::array<uchar, num_cols>& widths) {
        assert(widths[static_cast<size_t>(cols_traits::PRIMARY)] <= row_traits::PRIMARY_BITS);
        assert(widths[static_cast<size_t>(cols_traits::POINTER)] <= row_traits::POINTER_BITS);
        assert(widths[static_cast<size_t>(cols_traits::OFFSET)] <= row_traits::OFFSET_BITS);
        assert(widths[static_cast<size_t>(cols_traits::CHARACTER)] <= row_traits::CHARACTER_BITS);
    }

} __attribute__((packed));

// Column Sizes for rlbwt_columns supporting all ASCII characters
using rlbwt_columns_5 = rlbwt_row_bits<rlbwt_columns, 8, 8, 16, 8>; // (5 bytes)
using rlbwt_columns_6 = rlbwt_row_bits<rlbwt_columns, 8, 8, 24, 8>; // (6 bytes)
using rlbwt_columns_7 = rlbwt_row_bits<rlbwt_columns, 8, 12, 24, 12>; // (7 bytes)
using rlbwt_columns_8 = rlbwt_row_bits<rlbwt_columns, 8, 12, 32, 12>; // (8 bytes)
using rlbwt_columns_9 = rlbwt_row_bits<rlbwt_columns, 8, 16, 32, 16>; // (9 bytes)
using rlbwt_columns_10 = rlbwt_row_bits<rlbwt_columns, 8, 16, 40, 16>; // (10 bytes)
using rlbwt_columns_11 = rlbwt_row_bits<rlbwt_columns, 8, 20, 40, 20>; // (11 bytes)
using rlbwt_columns_12 = rlbwt_row_bits<rlbwt_columns, 8, 20, 48, 20>; // (12 bytes)
using rlbwt_columns_13 = rlbwt_row_bits<rlbwt_columns, 8, 24, 48, 24>; // (13 bytes)

// Column Sizes for rlbwt_columns_idx supporting all ASCII characters
using rlbwt_columns_idx_6 = rlbwt_row_bits<rlbwt_columns_idx, 8, 17, 15, 8>; // (6 bytes)
using rlbwt_columns_idx_7 = rlbwt_row_bits<rlbwt_columns_idx, 8, 21, 19, 8>; // (7 bytes)
using rlbwt_columns_idx_8 = rlbwt_row_bits<rlbwt_columns_idx, 8, 27, 25, 8>; // (8 bytes)
using rlbwt_columns_idx_9 = rlbwt_row_bits<rlbwt_columns_idx, 8, 29, 27, 12>; // (9 bytes)
using rlbwt_columns_idx_10 = rlbwt_row_bits<rlbwt_columns_idx, 8, 33, 31, 12>; // (10 bytes)
using rlbwt_columns_idx_11 = rlbwt_row_bits<rlbwt_columns_idx, 8, 35, 33, 12>; // (11 bytes)
using rlbwt_columns_idx_12 = rlbwt_row_bits<rlbwt_columns_idx, 8, 37, 35, 16>; // (12 bytes)
using rlbwt_columns_idx_13 = rlbwt_row_bits<rlbwt_columns_idx, 8, 41, 39, 16>; // (13 bytes)
using rlbwt_columns_idx_14 = rlbwt_row_bits<rlbwt_columns_idx, 8, 45, 43, 16>; // (14 bytes)
using rlbwt_columns_idx_15 = rlbwt_row_bits<rlbwt_columns_idx, 8, 47, 45, 20>; // (15 bytes)
using rlbwt_columns_idx_16 = rlbwt_row_bits<rlbwt_columns_idx, 8, 51, 49, 20>; // (16 bytes)

// Column Sizes for rlbwt_columns supporting DNA sequence characters
using dna_seq_columns_5 = rlbwt_row_bits<rlbwt_columns, 3, 10, 17, 10>; // (5 bytes)
using dna_seq_columns_6 = rlbwt_row_bits<rlbwt_columns, 3, 10, 25, 10>; // (6 bytes)
using dna_seq_columns_7 = rlbwt_row_bits<rlbwt_columns, 3, 14, 25, 14>; // (7 bytes)
using dna_seq_columns_8 = rlbwt_row_bits<rlbwt_columns, 3, 14, 33, 14>; // (8 bytes)
using dna_seq_columns_9 = rlbwt_row_bits<rlbwt_columns, 3, 18, 33, 18>; // (9 bytes)
using dna_seq_columns_10 = rlbwt_row_bits<rlbwt_columns, 3, 18, 41, 18>; // (10 bytes)
using dna_seq_columns_11 = rlbwt_row_bits<rlbwt_columns, 3, 22, 41, 22>; // (11 bytes)
using dna_seq_columns_12 = rlbwt_row_bits<rlbwt_columns, 3, 22, 49, 22>; // (12 bytes)
using dna_seq_columns_13 = rlbwt_row_bits<rlbwt_columns, 3, 26, 49, 26>; // (13 bytes)
using rlbwt_columns_default = dna_seq_columns_10;

// Column Sizes for rlbwt_columns_idx supporting all ASCII characters
using dna_seq_columns_idx_6 = rlbwt_row_bits<rlbwt_columns_idx, 3, 19, 17, 9>; // (6 bytes)
using dna_seq_columns_idx_7 = rlbwt_row_bits<rlbwt_columns_idx, 3, 23, 21, 9>; // (7 bytes)
using dna_seq_columns_idx_8 = rlbwt_row_bits<rlbwt_columns_idx, 3, 29, 27, 9>; // (8 bytes)
using dna_seq_columns_idx_9 = rlbwt_row_bits<rlbwt_columns_idx, 3, 31, 29, 13>; // (9 bytes)
using dna_seq_columns_idx_10 = rlbwt_row_bits<rlbwt_columns_idx, 3, 35, 33, 13>; // (10 bytes)
using dna_seq_columns_idx_11 = rlbwt_row_bits<rlbwt_columns_idx, 3, 37, 35, 13>; // (11 bytes)
using dna_seq_columns_idx_12 = rlbwt_row_bits<rlbwt_columns_idx, 3, 39, 37, 17>; // (12 bytes)
using dna_seq_columns_idx_13 = rlbwt_row_bits<rlbwt_columns_idx, 3, 43, 41, 17>; // (13 bytes)
using dna_seq_columns_idx_14 = rlbwt_row_bits<rlbwt_columns_idx, 3, 47, 45, 17>; // (14 bytes)
using dna_seq_columns_idx_15 = rlbwt_row_bits<rlbwt_columns_idx, 3, 49, 47, 21>; // (15 bytes)
using dna_seq_columns_idx_16 = rlbwt_row_bits<rlbwt_columns_idx, 3, 53, 51, 21>; // (16 bytes)
using rlbwt_columns_idx_default = dna_seq_columns_idx_13;

} // namespace orbit::rlbwt

namespace orbit {
template <>
struct move_row_traits<rlbwt::rlbwt_columns> : move_row_traits<rlbwt::rlbwt_columns_default> {};
template <>
struct move_row_traits<rlbwt::rlbwt_columns_idx> : move_row_traits<rlbwt::rlbwt_columns_idx_default> {};

template<>
struct table_row_for<rlbwt::rlbwt_columns> {
    using type = rlbwt::rlbwt_row<rlbwt::rlbwt_columns>;
};

template<>
struct table_row_for<rlbwt::rlbwt_columns_idx> {
    using type = rlbwt::rlbwt_row<rlbwt::rlbwt_columns_idx>;
};
} // namespace orbit

// ============================================= INVERTIBLE =============================================

namespace orbit {

template <typename columns_t, size_t C, size_t P, size_t PTR>
struct move_row_traits<rlbwt::invertible_rlbwt_row_bits<columns_t, C, P, PTR>> {
    static constexpr size_t CHARACTER_BITS = C;
    
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_FWD_BITS = PTR;
    static constexpr size_t POINTER_INV_BITS = PTR;
    
    static constexpr size_t FWD_INTERVAL_BITS = 1;
    static constexpr size_t INV_INTERVAL_BITS = 1;
};

} // namespace orbit

namespace orbit::rlbwt {

template <typename columns_t = rlbwt_columns>
struct invertible_rlbwt_row {
    // Sets num_cols, columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)
    using row_traits = move_row_traits<columns>;

    ulint primary : row_traits::PRIMARY_BITS;
    ulint pointer_fwd : row_traits::POINTER_FWD_BITS;
    ulint pointer_inv : row_traits::POINTER_INV_BITS;
    ulint fwd_interval : row_traits::FWD_INTERVAL_BITS;
    ulint inv_interval : row_traits::INV_INTERVAL_BITS;
    ulint character : row_traits::CHARACTER_BITS;

    invertible_rlbwt_row() = default;
    invertible_rlbwt_row(const std::array<ulint, num_cols>& values) {
        set(values);
    }

    template <columns col>
    void set(ulint val) {
        if constexpr (col == cols_traits::PRIMARY) primary = val;
        else if constexpr (col == cols_traits::POINTER_FWD) pointer_fwd = val;
        else if constexpr (col == cols_traits::POINTER_INV) pointer_inv = val;
        else if constexpr (col == cols_traits::FWD_INTERVAL) fwd_interval = val;
        else if constexpr (col == cols_traits::INV_INTERVAL) inv_interval = val;
        else if constexpr (col == cols_traits::CHARACTER) character = val;
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
        else if constexpr (col == cols_traits::CHARACTER) return character;
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
            static_cast<uchar>(row_traits::CHARACTER_BITS),
        };
        return widths;
    }

    static void assert_widths(const std::array<uchar, num_cols>& widths) {
        assert(widths[static_cast<size_t>(cols_traits::PRIMARY)] <= row_traits::PRIMARY_BITS);
        assert(widths[static_cast<size_t>(cols_traits::POINTER_FWD)] <= row_traits::POINTER_FWD_BITS);
        assert(widths[static_cast<size_t>(cols_traits::POINTER_INV)] <= row_traits::POINTER_INV_BITS);
        assert(widths[static_cast<size_t>(cols_traits::FWD_INTERVAL)] <= row_traits::FWD_INTERVAL_BITS);
        assert(widths[static_cast<size_t>(cols_traits::INV_INTERVAL)] <= row_traits::INV_INTERVAL_BITS);
        assert(widths[static_cast<size_t>(cols_traits::CHARACTER)] <= row_traits::CHARACTER_BITS);
    }

} __attribute__((packed));

// Column Sizes for rlbwt_columns supporting all ASCII characters
// Uses 2xPTR + 2 bits for fwd/inv interval
using invertible_rlbwt_cols_15 = invertible_rlbwt_row_bits<rlbwt_columns, 8, 16, 39>; // (13 bytes)
using invertible_rlbwt_cols_idx_17 = invertible_rlbwt_row_bits<rlbwt_columns_idx, 8, 42, 38>; // (16 bytes)
using invertible_dna_seq_columns_15 = invertible_rlbwt_row_bits<rlbwt_columns, 3, 19, 40>; // (13 bytes)
using invertible_dna_seq_cols_idx_17 = invertible_rlbwt_row_bits<rlbwt_columns_idx, 3, 45, 39>; // (16 bytes)

using invertible_rlbwt_columns_default = invertible_dna_seq_columns_15;
using invertible_rlbwt_columns_idx_default = invertible_dna_seq_cols_idx_17;

} // namespace orbit::rlbwt

namespace orbit {
template <>
struct move_row_traits<rlbwt::rlbwt_invertible_columns> : move_row_traits<rlbwt::invertible_rlbwt_columns_default> {};
template <>
struct move_row_traits<rlbwt::rlbwt_invertible_columns_idx> : move_row_traits<rlbwt::invertible_rlbwt_columns_idx_default> {};

template<>
struct table_row_for<rlbwt::rlbwt_invertible_columns> {
    using type = rlbwt::invertible_rlbwt_row<rlbwt::rlbwt_invertible_columns>;
};

template<>
struct table_row_for<rlbwt::rlbwt_invertible_columns_idx> {
    using type = rlbwt::invertible_rlbwt_row<rlbwt::rlbwt_invertible_columns_idx>;
};
} // namespace orbit

#endif /* end of include guard: _RLBWT_ROW_HPP */