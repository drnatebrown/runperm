#ifndef _MOVE_ROW_HPP
#define _MOVE_ROW_HPP

#include "common.hpp"
#include "move_columns.hpp"
#include <stdexcept>
#include <array>

// template <>
// struct MoveRowTraits<MoveCols> {
//     static constexpr size_t PRIMARY_BITS = BYTES_TO_BITS(LENGTH_BYTES);
//     static constexpr size_t POINTER_BITS = BYTES_TO_BITS(POINTER_BYTES);
//     static constexpr size_t OFFSET_BITS = BYTES_TO_BITS(OFFSET_BYTES);
// };

// template <>
// struct MoveRowTraits<MoveColsIdx> {
//     static constexpr size_t PRIMARY_BITS = BYTES_TO_BITS(START_BYTES);
//     static constexpr size_t POINTER_BITS = BYTES_TO_BITS(POINTER_BYTES);
//     static constexpr size_t OFFSET_BITS = BYTES_TO_BITS(OFFSET_BYTES);
// };

template <typename ColumnType>
struct MoveRowTraits;

template <typename Cols, size_t PRIMARY_BITS, size_t POINTER_BITS, size_t OFFSET_BITS>
struct MoveRowBits {};

template <typename Cols, size_t P, size_t PTR, size_t OFF>
struct MoveRowTraits<MoveRowBits<Cols, P, PTR, OFF>> {
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_BITS = PTR;
    static constexpr size_t OFFSET_BITS = OFF;
};

template <typename Columns = MoveCols>
struct MoveRow {
    using Columns = Columns;

    ulint primary : MoveRowTraits<Columns>::PRIMARY_BITS;
    ulint pointer : MoveRowTraits<Columns>::POINTER_BITS;
    ulint offset : MoveRowTraits<Columns>::OFFSET_BITS;

    MoveRow() = default;
    MoveRow(const std::array<ulint, Columns::NUM_COLS>& values) {
        set(values);
    }

    template <Columns Col>
    void set(ulint val) {
        if constexpr (Col == MoveColsTraits<Columns>::PRIMARY) primary = val;
        else if constexpr (Col == Columns::POINTER) pointer = val;
        else if constexpr (Col == Columns::OFFSET) offset = val;
    }
    template <size_t... Indices>
    void set(const std::array<ulint, Columns::NUM_COLS>& values, std::index_sequence<Indices...>) {
        (set<static_cast<Columns>(Indices)>(values[Indices]), ...);
    }
    void set(const std::array<ulint, Columns::NUM_COLS>& values) {
        set(values, std::make_index_sequence<Columns::NUM_COLS>{});
    }

    template <Columns Col>
    ulint get() const {
        if constexpr (Col == MoveColsTraits<Columns>::PRIMARY) return primary;
        else if constexpr (Col == Columns::POINTER) return pointer;
        else if constexpr (Col == Columns::OFFSET)   return offset;
        else throw std::invalid_argument("Invalid column");
    }
    template <size_t... Indices>
    std::array<ulint, Columns::NUM_COLS> get(std::index_sequence<Indices...>) const {
        return {get<static_cast<Columns>(Indices)>()...};
    }
    std::array<ulint, Columns::NUM_COLS> get() const {
        return get(std::make_index_sequence<Columns::NUM_COLS>{});
    }
} __attribute__((packed));

// Column Sizes for MoveCols
using MoveCol_4 = MoveRowBits<MoveCols, 8, 16, 8>;     // (4 bytes)
using MoveCol_5 = MoveRowBits<MoveCols, 8, 24, 8>;     // (5 bytes)
using MoveCol_6 = MoveRowBits<MoveCols, 12, 24, 12>; // (6 bytes)
using MoveCol_7 = MoveRowBits<MoveCols, 12, 32, 12>; // (7 bytes)
using MoveCol_8 = MoveRowBits<MoveCols, 16, 32, 16>; // (8 bytes)
using MoveCol_9 = MoveRowBits<MoveCols, 16, 40, 16>; // (9 bytes)
using MoveCol_10 = MoveRowBits<MoveCols, 20, 40, 20>; // (10 bytes)
using MoveCol_11 = MoveRowBits<MoveCols, 20, 48, 20>; // (11 bytes)
using MoveCol_12 = MoveRowBits<MoveCols, 24, 48, 24>; // (12 bytes)
using MoveColDefault = MoveCol_9;

// Column Sizes for MoveRowIdx
using MoveColIdx_5 = MoveRowBits<MoveColsIdx, 17, 15, 8>; // (5 bytes)
using MoveColIdx_6 = MoveRowBits<MoveColsIdx, 21, 19, 8>; // (6 bytes)
using MoveColIdx_7 = MoveRowBits<MoveColsIdx, 27, 25, 8>; // (7 bytes)
using MoveColIdx_8 = MoveRowBits<MoveColsIdx, 29, 27, 12>; // (8 bytes)
using MoveColIdx_9 = MoveRowBits<MoveColsIdx, 33, 31, 12>; // (9 bytes)
using MoveColIdx_10 = MoveRowBits<MoveColsIdx, 35, 33, 12>; // (10 bytes)
using MoveColIdx_11 = MoveRowBits<MoveColsIdx, 37, 35, 16>; // (11 bytes)
using MoveColIdx_12 = MoveRowBits<MoveColsIdx, 41, 39, 16>; // (12 bytes)
using MoveColIdx_13 = MoveRowBits<MoveColsIdx, 45, 43, 16>; // (13 bytes)
using MoveColIdx_14 = MoveRowBits<MoveColsIdx, 47, 45, 20>; // (14 bytes)
using MoveColIdx_15 = MoveRowBits<MoveColsIdx, 51, 49, 20>; // (15 bytes)
using MoveColIdxDefault = MoveColIdx_12;

template <>
struct MoveRowTraits<MoveCols> : MoveRowTraits<MoveColDefault> {};

template <>
struct MoveRowTraits<MoveColsIdx> : MoveRowTraits<MoveColIdxDefault> {};

#endif /* end of include guard: _MOVE_ROW_HPP */