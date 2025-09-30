#ifndef _MOVE_ROW_HPP
#define _MOVE_ROW_HPP

#include "common.hpp"
#include "move/move_columns.hpp"
#include <stdexcept>
#include <array>

template <typename ColumnsType>
struct MoveRowTraits;

template <typename ColsType, size_t PRIMARY_BITS, size_t POINTER_BITS, size_t OFFSET_BITS>
struct MoveRowBits {};

template <typename ColsType, size_t P, size_t PTR, size_t OFF>
struct MoveRowTraits<MoveRowBits<ColsType, P, PTR, OFF>> {
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_BITS = PTR;
    static constexpr size_t OFFSET_BITS = OFF;
};

template <typename ColumnsType = MoveCols>
struct MoveRow {
    MOVE_CLASS_TRAITS(ColumnsType)
    using RowTraits = MoveRowTraits<Columns>;

    ulint primary : RowTraits::PRIMARY_BITS;
    ulint pointer : RowTraits::POINTER_BITS;
    ulint offset : RowTraits::OFFSET_BITS;

    MoveRow() = default;
    MoveRow(const std::array<ulint, NUM_COLS>& values) {
        set(values);
    }

    template <Columns Col>
    void set(ulint val) {
        if constexpr (Col == ColsTraits::PRIMARY) primary = val;
        else if constexpr (Col == ColsTraits::POINTER) pointer = val;
        else if constexpr (Col == ColsTraits::OFFSET) offset = val;
    }
    template <size_t... Indices>
    void set(const std::array<ulint, NUM_COLS>& values, std::index_sequence<Indices...>) {
        (set<static_cast<Columns>(Indices)>(values[Indices]), ...);
    }
    void set(const std::array<ulint, NUM_COLS>& values) {
        set(values, std::make_index_sequence<NUM_COLS>{});
    }

    template <Columns Col>
    ulint get() const {
        if constexpr (Col == ColsTraits::PRIMARY) return primary;
        else if constexpr (Col == ColsTraits::POINTER) return pointer;
        else if constexpr (Col == ColsTraits::OFFSET)   return offset;
        else throw std::invalid_argument("Invalid column");
    }
    template <size_t... Indices>
    std::array<ulint, NUM_COLS> get(std::index_sequence<Indices...>) const {
        return {get<static_cast<Columns>(Indices)>()...};
    }
    std::array<ulint, NUM_COLS> get() const {
        return get(std::make_index_sequence<NUM_COLS>{});
    }

    static void assert_widths(const std::array<uchar, NUM_COLS>& widths) {
        assert(widths[static_cast<size_t>(ColsTraits::PRIMARY)] <= RowTraits::PRIMARY_BITS);
        assert(widths[static_cast<size_t>(ColsTraits::POINTER)] <= RowTraits::POINTER_BITS);
        assert(widths[static_cast<size_t>(ColsTraits::OFFSET)] <= RowTraits::OFFSET_BITS);
    }

} __attribute__((packed));

// Column Sizes for MoveCols
using MoveCols_4 = MoveRowBits<MoveCols, 8, 16, 8>;     // (4 bytes)
using MoveCols_5 = MoveRowBits<MoveCols, 8, 24, 8>;     // (5 bytes)
using MoveCols_6 = MoveRowBits<MoveCols, 12, 24, 12>; // (6 bytes)
using MoveCols_7 = MoveRowBits<MoveCols, 12, 32, 12>; // (7 bytes)
using MoveCols_8 = MoveRowBits<MoveCols, 16, 32, 16>; // (8 bytes)
using MoveCols_9 = MoveRowBits<MoveCols, 16, 40, 16>; // (9 bytes)
using MoveCols_10 = MoveRowBits<MoveCols, 20, 40, 20>; // (10 bytes)
using MoveCols_11 = MoveRowBits<MoveCols, 20, 48, 20>; // (11 bytes)
using MoveCols_12 = MoveRowBits<MoveCols, 24, 48, 24>; // (12 bytes)
using MoveColsDefault = MoveCols_9;

// Column Sizes for MoveRowIdx
using MoveColsIdx_5 = MoveRowBits<MoveColsIdx, 17, 15, 8>; // (5 bytes)
using MoveColsIdx_6 = MoveRowBits<MoveColsIdx, 21, 19, 8>; // (6 bytes)
using MoveColsIdx_7 = MoveRowBits<MoveColsIdx, 27, 25, 8>; // (7 bytes)
using MoveColsIdx_8 = MoveRowBits<MoveColsIdx, 29, 27, 12>; // (8 bytes)
using MoveColsIdx_9 = MoveRowBits<MoveColsIdx, 33, 31, 12>; // (9 bytes)
using MoveColsIdx_10 = MoveRowBits<MoveColsIdx, 35, 33, 12>; // (10 bytes)
using MoveColsIdx_11 = MoveRowBits<MoveColsIdx, 37, 35, 16>; // (11 bytes)
using MoveColsIdx_12 = MoveRowBits<MoveColsIdx, 41, 39, 16>; // (12 bytes)
using MoveColsIdx_13 = MoveRowBits<MoveColsIdx, 45, 43, 16>; // (13 bytes)
using MoveColsIdx_14 = MoveRowBits<MoveColsIdx, 47, 45, 20>; // (14 bytes)
using MoveColsIdx_15 = MoveRowBits<MoveColsIdx, 51, 49, 20>; // (15 bytes)
using MoveColsIdxDefault = MoveColsIdx_12;

template <>
struct MoveRowTraits<MoveCols> : MoveRowTraits<MoveColsDefault> {};

template <>
struct MoveRowTraits<MoveColsIdx> : MoveRowTraits<MoveColsIdxDefault> {};

#endif /* end of include guard: _MOVE_ROW_HPP */