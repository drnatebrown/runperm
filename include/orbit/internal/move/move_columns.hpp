#ifndef _MOVE_COLUMNS_HPP
#define _MOVE_COLUMNS_HPP

#include "orbit/common.hpp"

namespace orbit {

/* ALL MOVE STRUCTURE MUST REPRESENT AT LEAST ONE OF THESE COLUMNS DEFS */

/* Stores permutation π(1..n) in O(r)-space, where r is the number of positions such that
either i == 0 or π(i-1) != π(i) - 1. An interval of length l is any range where π(i) is stored, and
otherwise π(j), for i <=j < (i - l), is computed as π(j) = π(i) + (j - i). Here, we use interval/offset notation. */
enum class move_columns { 
    LENGTH, // Length of move interval
    POINTER, // Where i is start of current interval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    COUNT // Helper to get the number of columns
};

/* Above, but uses absolute indices instead of relative intervals */
enum class move_columns_idx { 
    START, // i is the start of the interval
    POINTER, // Where i is start of current inverval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    COUNT // Helper to get the number of columns
};

template<bool is_relative>
struct move_position;

template <>
struct move_position<true> {
    struct type {
        ulint interval = 0;
        ulint offset = 0;

        bool operator==(const type& other) const {
            return interval == other.interval && offset == other.offset;
        }
        bool operator!=(const type& other) const {
            return !(*this == other);
        }
    };
};

template <>
struct move_position<false> {
    struct type {
        ulint interval = 0;
        ulint offset = 0;
        ulint idx = 0;

        bool operator==(const type& other) const {
            return interval == other.interval && offset == other.offset;
        }
        bool operator!=(const type& other) const {
            return !(*this == other);
        }
    };
};

template <typename columns_t>
struct move_cols_traits;
/* Needs to define:
   - RELATIVE: bool // whether position is relative or absolute
   - PRIMARY: enum class Columns
   - POINTER: enum class Columns
   - OFFSET: enum class Columns
   - NUM_COLS: size_t
   - position: struct with interval, offset, and idx (if ABSOLUTE)
*/

template <>
struct move_cols_traits<move_columns> {
    static constexpr bool RELATIVE = true;
    static constexpr bool INVERTIBLE = false;
    static constexpr move_columns PRIMARY = move_columns::LENGTH;

    static constexpr move_columns LENGTH = move_columns::LENGTH;
    static constexpr move_columns POINTER = move_columns::POINTER;
    static constexpr move_columns OFFSET = move_columns::OFFSET;
    static constexpr size_t NUM_COLS = static_cast<size_t>(move_columns::COUNT);

    using position = move_position<RELATIVE>::type;
};

template <>
struct move_cols_traits<move_columns_idx> {
    static constexpr bool RELATIVE = false;
    static constexpr bool INVERTIBLE = false;
    static constexpr move_columns_idx PRIMARY = move_columns_idx::START;

    static constexpr move_columns_idx START = move_columns_idx::START;
    static constexpr move_columns_idx POINTER = move_columns_idx::POINTER;
    static constexpr move_columns_idx OFFSET = move_columns_idx::OFFSET;
    static constexpr size_t NUM_COLS = static_cast<size_t>(move_columns_idx::COUNT);
    
    using position = move_position<RELATIVE>::type;
};

/* Allows for switching between Length and Start columns, if passed one type but want relative or absolute positions.
   Specialized move columns should also implement this
*/
template<typename base_columns>
struct column_switcher;

template<>
struct column_switcher<move_columns> {
    using relative = move_columns;
    using absolute = move_columns_idx;
};

template<>
struct column_switcher<move_columns_idx> {
    using relative = move_columns;
    using absolute = move_columns_idx;
};

// Helper alias for easy switching
template<typename base_columns, bool use_absolute_positions>
using switch_columns = std::conditional_t<use_absolute_positions,
    typename column_switcher<base_columns>::absolute,
    typename column_switcher<base_columns>::relative>;

/* Code below is only useful to interact with permutation, so we can have base column oblivious extended types */
// ADL detector: is C an extended enum (has friend columns_parent(C))?
template<class c, class = void>
struct is_columns_extended : std::false_type {};

template<class c>
struct is_columns_extended<c, std::void_t<decltype(columns_parent(std::declval<c>()))>>
    : std::true_type {};

// Get the parent RunDataColumns<RunData,Base> for extended enums via ADL
template<class c>
using columns_parent_t = std::remove_pointer_t<decltype(columns_parent(std::declval<c>()))>;

// Resolver can now resolve extended columns types, or just use existing specializations
template<class c, bool = is_columns_extended<c>::value>
struct resolve_cols_traits { using type = move_cols_traits<c>; };

// ================================ INVERTIBLE MOVE STRUCTURES ================================
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
struct move_cols_traits<invertible_columns> {
    static constexpr bool RELATIVE = true;
    static constexpr bool INVERTIBLE = true;
    static constexpr invertible_columns PRIMARY = invertible_columns::LENGTH;

    static constexpr invertible_columns LENGTH = invertible_columns::LENGTH;
    static constexpr invertible_columns POINTER_FWD = invertible_columns::POINTER_FWD;
    static constexpr invertible_columns POINTER_INV = invertible_columns::POINTER_INV;
    static constexpr invertible_columns FWD_INTERVAL = invertible_columns::FWD_INTERVAL;
    static constexpr invertible_columns INV_INTERVAL = invertible_columns::INV_INTERVAL;
    static constexpr size_t NUM_COLS = static_cast<size_t>(invertible_columns::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct move_cols_traits<invertible_columns_idx> {
    static constexpr bool RELATIVE = false;
    static constexpr bool INVERTIBLE = true;
    static constexpr invertible_columns_idx PRIMARY = invertible_columns_idx::START;

    static constexpr invertible_columns_idx START = invertible_columns_idx::START;
    static constexpr invertible_columns_idx POINTER_FWD = invertible_columns_idx::POINTER_FWD;
    static constexpr invertible_columns_idx POINTER_INV = invertible_columns_idx::POINTER_INV;
    static constexpr invertible_columns_idx FWD_INTERVAL = invertible_columns_idx::FWD_INTERVAL;
    static constexpr invertible_columns_idx INV_INTERVAL = invertible_columns_idx::INV_INTERVAL;
    static constexpr size_t NUM_COLS = static_cast<size_t>(invertible_columns_idx::COUNT);

    using position = move_position<RELATIVE>::type;
};

template<>
struct column_switcher<invertible_columns> {
    using relative = invertible_columns;
    using absolute = invertible_columns_idx;
};

template<>
struct column_switcher<invertible_columns_idx> {
    using relative = invertible_columns;
    using absolute = invertible_columns_idx;
};
} // namespace orbit

#endif /* end of include guard: _MOVE_COLUMNS_HPP */