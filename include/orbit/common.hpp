#ifndef _COMMON_HPP
#define _COMMON_HPP
#include <stddef.h>
#include <utility>
#include <vector>
#include <array>
#include <cmath>
#include <fstream>
#include <tuple>

namespace orbit {

typedef unsigned char uchar;
typedef unsigned long int ulint;

// LEAVE UNCHANGED
inline constexpr size_t MAX_ALPHABET_SIZE = 256;
inline constexpr size_t RW_BYTES = 5;

// CONFIGURABLE ===============================================================
inline constexpr uchar TERMINATOR = 0;
inline constexpr uchar SEPARATOR = 1;

inline constexpr size_t START_BYTES = 5;
inline constexpr size_t LENGTH_BYTES = 2;
inline constexpr size_t POINTER_BYTES = 4;
inline constexpr size_t OFFSET_BYTES = 2;
inline constexpr size_t CHARACTER_BYTES = 1;
inline constexpr size_t DEFAULT_BYTES = 4;
// END CONFIGURABLE ===========================================================

inline constexpr size_t VERSION_MAJOR = 1;
inline constexpr size_t VERSION_MINOR = 0;
inline constexpr size_t VERSION_PATCH = 0;

inline size_t serialize_version(std::ostream& out) {
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

inline std::tuple<size_t, size_t, size_t> load_version(std::istream& is) {
    size_t major, minor, patch;
    is.read(reinterpret_cast<char *>(&major), sizeof(major));
    is.read(reinterpret_cast<char *>(&minor), sizeof(minor));
    is.read(reinterpret_cast<char *>(&patch), sizeof(patch));
    return {major, minor, patch};
}

inline constexpr size_t bytes_to_bits(size_t bytes) { return bytes * 8; }
inline constexpr size_t bits_to_bytes(size_t bits) { return bits / 8; }
#define num_bits_type(type) bytes_to_bits(sizeof(type))
inline constexpr size_t ceil_div(size_t num, size_t den) { return (num + den - 1) / den; }
inline constexpr size_t pow2(size_t bits) { return 1ULL << bits; }
inline constexpr size_t max_val(size_t bits) { return pow2(bits) - 1; }
inline constexpr size_t mask(size_t bits) { return max_val(bits); }

inline constexpr const char MOVE_STRUCTURE_EXTENSION[] = ".move";

inline constexpr uchar bit_width(ulint value) {
    return value == 0 ? 1 : 64 - __builtin_clzll(value);
}

// ENUM HELPERS ===============================================================
// ENUM REPRESENTS COLUMNS, USE ENUM HELPERS TO ENFORCE STRUCTURE
// in common.hpp (or an internal traits header
template<class, class = void>
struct has_count_enumerator : std::false_type {};

template<class E>
struct has_count_enumerator<E, std::void_t<decltype(E::COUNT)>> : std::true_type {};

template<class E>
inline constexpr size_t to_index(E e) noexcept { return static_cast<size_t>(e); }

template<class E>
inline constexpr size_t num_columns() noexcept {
    return static_cast<size_t>(E::COUNT);
}

template<class E, typename T = ulint>
using columns_tuple = std::array<T, num_columns<E>()>;

#define MOVE_CLASS_TRAITS(columns_param) \
    using columns = columns_param; \
    template<typename C> \
    using cols_traits_for = typename resolve_cols_traits<C>::type;  \
    using cols_traits = typename resolve_cols_traits<columns>::type;  \
    static constexpr size_t num_cols = cols_traits::NUM_COLS; \
    template<typename E> static constexpr columns to_cols(E e) { return static_cast<columns>(e); }

} // namespace orbit

// Defines a macro to generate an enum class named <enum_name> with specified fields, 
// and appends COUNT as the last enumerator for sizing.
// Usage: DEFINE_ORBIT_COLUMNS(MyEnum, FIELD1, FIELD2, FIELD3)
#define DEFINE_ORBIT_COLUMNS(enum_name, ...) \
    enum class enum_name { __VA_ARGS__, COUNT };

#endif /* end of include guard: _COMMON_HPP */