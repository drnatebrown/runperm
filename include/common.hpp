#ifndef _COMMON_HPP
#define _COMMON_HPP

// LEAVE UNCHANGED
#define ALPHABET_SIZE 256
#define RW_BYTES 5

// CONFIGURABLE
#define TERMINATOR 1
#define SEPARATOR 2

#define START_BYTES 5
#define LENGTH_BYTES 2
#define POINTER_BYTES 4
#define OFFSET_BYTES 2

#define BYTES_TO_BITS(bytes) ((bytes) * 8)
#define BITS_TO_BYTES(bits) ((bits) / 8)
#define NUM_BITS(type) (BYTES_TO_BITS(sizeof(type))) 
#define CEIL_DIV(num, den) ((num + den - 1) / den)

#define MOVE_STRUCTURE_EXTENSION ".move"

#define MOVE_CLASS_TRAITS(ColumnsParam) \
    using Columns = ColumnsParam; \
    using ColsTraits = MoveColsTraits<Columns>; \
    static constexpr size_t NUM_COLS = ColsTraits::NUM_COLS;

typedef unsigned char uchar;
typedef unsigned long int ulint;

constexpr uchar bit_width(ulint value) {
    return value == 0 ? 1 : 64 - __builtin_clzll(value);
}
#endif /* end of include guard: _COMMON_HPP */