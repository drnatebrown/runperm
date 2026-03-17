#ifndef _COMMON_HPP
#define _COMMON_HPP
#include <stddef.h>
#include <utility>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <tuple>

#define VERSION_MAJOR 0
#define VERSION_MINOR 4
#define VERSION_PATCH 0
#define VERSION_STRING VERSION_MAJOR "." VERSION_MINOR "." VERSION_PATCH

size_t serialize_version(std::ostream& out) {
    size_t written_bytes = 0;
    size_t major = VERSION_MAJOR;
    size_t minor = VERSION_MINOR;
    size_t patch = VERSION_PATCH;
    out.write(reinterpret_cast<char *>(&major), sizeof(major));
    written_bytes += sizeof(major);
    out.write(reinterpret_cast<char *>(&minor), sizeof(minor));
    written_bytes += sizeof(minor);
    out.write(reinterpret_cast<char *>(&patch), sizeof(patch));
    written_bytes += sizeof(patch);
    return written_bytes;
}

std::tuple<size_t, size_t, size_t> load_version(std::istream& is) {
    size_t major, minor, patch;
    is.read(reinterpret_cast<char *>(&major), sizeof(major));
    is.read(reinterpret_cast<char *>(&minor), sizeof(minor));
    is.read(reinterpret_cast<char *>(&patch), sizeof(patch));
    return {major, minor, patch};
}

// LEAVE UNCHANGED
#define MAX_ALPHABET_SIZE 256
#define RW_BYTES 5

// CONFIGURABLE
#define TERMINATOR 0
#define SEPARATOR 1

#define START_BYTES 5
#define LENGTH_BYTES 2
#define POINTER_BYTES 4
#define OFFSET_BYTES 2
#define CHARACTER_BYTES 1
#define DEFAULT_BYTES 4

#define BYTES_TO_BITS(bytes) ((bytes) * 8)
#define BITS_TO_BYTES(bits) ((bits) / 8)
#define NUM_BITS(type) (BYTES_TO_BITS(sizeof(type))) 
#define CEIL_DIV(num, den) ((num + den - 1) / den)
#define POW2(bits) (1ULL << (bits))
#define MAX_VAL(bits) (POW2(bits) - 1)
#define MASK(bits) MAX_VAL(bits)

#define MOVE_STRUCTURE_EXTENSION ".move"

typedef unsigned char uchar;
typedef unsigned long int ulint;

constexpr uchar bit_width(ulint value) {
    return value == 0 ? 1 : 64 - __builtin_clzll(value);
}

// ENUM REPRESENTS COLUMNS, USE ENUM HELPERS TO ENFORCE STRUCTURE
// in common.hpp (or an internal traits header)
template<class, class = void>
struct has_count_enumerator : std::false_type {};

template<class E>
struct has_count_enumerator<E, std::void_t<decltype(E::COUNT)>> : std::true_type {};

template<class E>
constexpr size_t to_index(E e) noexcept { return static_cast<size_t>(e); }

template<class E>
constexpr size_t num_columns() noexcept {
    return static_cast<size_t>(E::COUNT);
}

template<class E, typename T = ulint>
using ColumnsTuple = std::array<T, num_columns<E>()>;

// Defines a macro to generate an enum class named <enum_name> with specified fields, 
// and appends COUNT as the last enumerator for sizing.
// Usage: DEFINE_COLUMNS(MyEnum, FIELD1, FIELD2, FIELD3)
#define DEFINE_COLUMNS(enum_name, ...) \
    enum class enum_name { __VA_ARGS__, COUNT };

#define MOVE_CLASS_TRAITS(ColumnsParam) \
    using Columns = ColumnsParam; \
    template<typename C> \
    using ColsTraitsFor = typename ResolveColsTraits<C>::type;  \
    using ColsTraits = typename ResolveColsTraits<Columns>::type;  \
    static constexpr size_t NumCols = ColsTraits::NUM_COLS; \
    template<typename E> static constexpr Columns to_cols(E e) { return static_cast<Columns>(e); }

inline std::pair<std::vector<uchar>, std::vector<ulint>> bwt_to_rlbwt(const std::vector<uchar> &bwt_chars) {
    std::vector<uchar> rlbwt_chars;
    std::vector<ulint> rlbwt_run_lengths;
    for (size_t i = 0; i < bwt_chars.size(); ++i) {
        if (i == 0 || bwt_chars[i] != bwt_chars[i - 1]) {
            rlbwt_chars.push_back(bwt_chars[i]);
            rlbwt_run_lengths.push_back(1);
        } else {
            ++rlbwt_run_lengths.back();
        }
    }
    return {rlbwt_chars, rlbwt_run_lengths};
}

#endif /* end of include guard: _COMMON_HPP */