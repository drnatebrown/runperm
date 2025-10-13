#ifndef _RLBWT_ROW_HPP
#define _RLBWT_ROW_HPP

#include "internal/common.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"
#include "internal/move/move_row.hpp"
#include <stdexcept>
#include <array>

// RLBWTRowBits takes bit values directly
template <size_t CHARACTER_BITS, size_t PRIMARY_BITS, size_t POINTER_BITS, size_t OFFSET_BITS>
struct RLBWTRowBits {
};

// Specialization for MoveRowTraits
template <size_t C, size_t P, size_t PTR, size_t OFF>
struct MoveRowTraits<RLBWTRowBits<C, P, PTR, OFF>> {
    static constexpr size_t CHARACTER_BITS = C;
    static constexpr size_t PRIMARY_BITS = P;
    static constexpr size_t POINTER_BITS = PTR;
    static constexpr size_t OFFSET_BITS = OFF;
};

template <typename ColumnsType = RLBWTCols>
struct RLBWTRow {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(ColumnsType)
    using RowTraits = MoveRowTraits<Columns>;

    ulint primary : RowTraits::PRIMARY_BITS;
    ulint pointer : RowTraits::POINTER_BITS;
    ulint offset : RowTraits::OFFSET_BITS;
    ulint character : RowTraits::CHARACTER_BITS;

    RLBWTRow() = default;
    RLBWTRow(const std::array<ulint, NumCols>& values) {
        set(values);
    }

    template <Columns Col>
    void set(ulint val) {
        if constexpr (Col == ColsTraits::PRIMARY) primary = val;
        else if constexpr (Col == ColsTraits::POINTER) pointer = val;
        else if constexpr (Col == ColsTraits::OFFSET) offset = val;
        else if constexpr (Col == ColsTraits::CHARACTER) character = val;
    }
    template <size_t... Indices>
    void set(const std::array<ulint, NumCols>& values, std::index_sequence<Indices...>) {
        (set<static_cast<Columns>(Indices)>(values[Indices]), ...);
    }
    void set(const std::array<ulint, NumCols>& values) {
        set(values, std::make_index_sequence<NumCols>{});
    }

    template <Columns Col>
    ulint get() const {
        if constexpr (Col == ColsTraits::PRIMARY) return primary;
        else if constexpr (Col == ColsTraits::POINTER) return pointer;
        else if constexpr (Col == ColsTraits::OFFSET)   return offset;
        else if constexpr (Col == ColsTraits::CHARACTER) return character;
        else throw std::invalid_argument("Invalid column");
    }
    template <size_t... Indices>
    std::array<ulint, NumCols> get(std::index_sequence<Indices...>) const {
        return {get<static_cast<Columns>(Indices)>()...};
    }
    std::array<ulint, NumCols> get() const {
        return get(std::make_index_sequence<NumCols>{});
    }

    static void assert_widths(const std::array<uchar, NumCols>& widths) {
        assert(widths[static_cast<size_t>(ColsTraits::PRIMARY)] <= RowTraits::PRIMARY_BITS);
        assert(widths[static_cast<size_t>(ColsTraits::POINTER)] <= RowTraits::POINTER_BITS);
        assert(widths[static_cast<size_t>(ColsTraits::OFFSET)] <= RowTraits::OFFSET_BITS);
        assert(widths[static_cast<size_t>(ColsTraits::CHARACTER)] <= RowTraits::CHARACTER_BITS);
    }

} __attribute__((packed));

// Column Sizes for RLBWTCols supporting all ASCII characters
using RLBWTCols_5 = RLBWTRowBits<8, 8, 16, 8>; // (5 bytes)
using RLBWTCols_6 = RLBWTRowBits<8, 8, 24, 8>; // (6 bytes)
using RLBWTCols_7 = RLBWTRowBits<8, 12, 24, 12>; // (7 bytes)
using RLBWTCols_8 = RLBWTRowBits<8, 12, 32, 12>; // (8 bytes)
using RLBWTCols_9 = RLBWTRowBits<8, 16, 32, 16>; // (9 bytes)
using RLBWTCols_10 = RLBWTRowBits<8, 16, 40, 16>; // (10 bytes)
using RLBWTCols_11 = RLBWTRowBits<8, 20, 40, 20>; // (11 bytes)
using RLBWTCols_12 = RLBWTRowBits<8, 20, 48, 20>; // (12 bytes)
using RLBWTCols_13 = RLBWTRowBits<8, 24, 48, 24>; // (13 bytes)

// Column Sizes for RLBWTColsIdx supporting all ASCII characters
using RLBWTColsIdx_6 = RLBWTRowBits<8, 17, 15, 8>; // (6 bytes)
using RLBWTColsIdx_7 = RLBWTRowBits<8, 21, 19, 8>; // (7 bytes)
using RLBWTColsIdx_8 = RLBWTRowBits<8, 27, 25, 8>; // (8 bytes)
using RLBWTColsIdx_9 = RLBWTRowBits<8, 29, 27, 12>; // (9 bytes)
using RLBWTColsIdx_10 = RLBWTRowBits<8, 33, 31, 12>; // (10 bytes)
using RLBWTColsIdx_11 = RLBWTRowBits<8, 35, 33, 12>; // (11 bytes)
using RLBWTColsIdx_12 = RLBWTRowBits<8, 37, 35, 16>; // (12 bytes)
using RLBWTColsIdx_13 = RLBWTRowBits<8, 41, 39, 16>; // (13 bytes)
using RLBWTColsIdx_14 = RLBWTRowBits<8, 45, 43, 16>; // (14 bytes)
using RLBWTColsIdx_15 = RLBWTRowBits<8, 47, 45, 20>; // (15 bytes)
using RLBWTColsIdx_16 = RLBWTRowBits<8, 51, 49, 20>; // (16 bytes)

// Column Sizes for RLBWTCols supporting DNA sequence characters
using DNASeqCols_5 = RLBWTRowBits<3, 10, 17, 10>; // (5 bytes)
using DNASeqCols_6 = RLBWTRowBits<3, 10, 25, 10>; // (6 bytes)
using DNASeqCols_7 = RLBWTRowBits<3, 14, 25, 14>; // (7 bytes)
using DNASeqCols_8 = RLBWTRowBits<3, 14, 33, 14>; // (8 bytes)
using DNASeqCols_9 = RLBWTRowBits<3, 18, 33, 18>; // (9 bytes)
using DNASeqCols_10 = RLBWTRowBits<3, 18, 41, 18>; // (10 bytes)
using DNASeqCols_11 = RLBWTRowBits<3, 22, 41, 22>; // (11 bytes)
using DNASeqCols_12 = RLBWTRowBits<3, 22, 49, 22>; // (12 bytes)
using DNASeqCols_13 = RLBWTRowBits<3, 26, 49, 26>; // (13 bytes)
using RLBWTColsDefault = DNASeqCols_10;

// Column Sizes for RLBWTColsIdx supporting all ASCII characters
using DNASeqColsIdx6 = RLBWTRowBits<3, 19, 17, 9>; // (6 bytes)
using DNASeqColsIdx7 = RLBWTRowBits<3, 23, 21, 9>; // (7 bytes)
using DNASeqColsIdx8 = RLBWTRowBits<3, 29, 27, 9>; // (8 bytes)
using DNASeqColsIdx9 = RLBWTRowBits<3, 31, 29, 13>; // (9 bytes)
using DNASeqColsIdx10 = RLBWTRowBits<3, 35, 33, 13>; // (10 bytes)
using DNASeqColsIdx11 = RLBWTRowBits<3, 37, 35, 13>; // (11 bytes)
using DNASeqColsIdx12 = RLBWTRowBits<3, 39, 37, 17>; // (12 bytes)
using DNASeqColsIdx13 = RLBWTRowBits<3, 43, 41, 17>; // (13 bytes)
using DNASeqColsIdx14 = RLBWTRowBits<3, 47, 45, 17>; // (14 bytes)
using DNASeqColsIdx15 = RLBWTRowBits<3, 49, 47, 21>; // (15 bytes)
using DNASeqColsIdx16 = RLBWTRowBits<3, 53, 51, 21>; // (16 bytes)
using RLBWTColsIdxDefault = DNASeqColsIdx13;

template <>
struct MoveRowTraits<RLBWTCols> : MoveRowTraits<RLBWTColsDefault> {};
template <>
struct MoveRowTraits<RLBWTColsIdx> : MoveRowTraits<RLBWTColsIdxDefault> {};

#endif /* end of include guard: _MOVE_ROW_HPP */