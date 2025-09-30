#ifndef _MOVE_COLUMNS_HPP
#define _MOVE_COLUMNS_HPP

#include "common.hpp"

/* ALL MOVE STRUCTURE MUST REPRESENT AT LEAST ONE OF THESE COLUMNS DEFS */

/* Stores permutation π(1..n) in O(r)-space, where r is the number of positions such that
either i == 0 or π(i-1) != π(i) - 1. An interval of length l is any range where π(i) is stored, and
otherwise π(j), for i <=j < (i - l), is computed as π(j) = π(i) + (j - i). Here, we use interval/offset notation. */
enum class MoveCols { 
    LENGTH, // Length of move interval
    POINTER, // Where i is start of current interval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    NUM_COLS // Helper to get the number of columns
};

/* Above, but uses absolute indices instead of relative intervals */
enum class MoveColsIdx { 
    START, // i is the start of the interval
    POINTER, // Where i is start of current inverval, the interval containing π(i)
    OFFSET, // Offset of π(i) within that interval
    NUM_COLS // Helper to get the number of columns
};

template <typename ColumnsType>
struct MoveColsTraits;
/* Needs to define:
   - ABSOLUTE: bool
   - RELATIVE: bool
   - PRIMARY: enum class Columns
   - POINTER: enum class Columns
   - OFFSET: enum class Columns
   - NUM_COLS: size_t
   - Position: struct with interval, offset, and idx (if ABSOLUTE)
   // TODO couple position type to absolute/relative
*/

template <>
struct MoveColsTraits<MoveCols> {
    static constexpr bool ABSOLUTE = false;
    static constexpr bool RELATIVE = true;
    static constexpr MoveCols PRIMARY = MoveCols::LENGTH;
    static constexpr MoveCols POINTER = MoveCols::POINTER;
    static constexpr MoveCols OFFSET = MoveCols::OFFSET;
    static constexpr size_t NUM_COLS = static_cast<size_t>(MoveCols::NUM_COLS);

    struct Position {
        ulint interval = 0;
        ulint offset = 0;
    };

    static constexpr bool HAS_LENGTH = true;
    static constexpr bool HAS_START = false;
    static constexpr MoveCols LENGTH = MoveCols::LENGTH;
};

template <>
struct MoveColsTraits<MoveColsIdx> {
    static constexpr bool ABSOLUTE = true;
    static constexpr bool RELATIVE = false;
    static constexpr MoveColsIdx PRIMARY = MoveColsIdx::START;
    static constexpr MoveColsIdx POINTER = MoveColsIdx::POINTER;
    static constexpr MoveColsIdx OFFSET = MoveColsIdx::OFFSET;
    static constexpr size_t NUM_COLS = static_cast<size_t>(MoveColsIdx::NUM_COLS);
    
    struct Position {
        ulint interval = 0;
        ulint offset = 0;
        ulint idx = 0;
    };

    static constexpr bool HAS_LENGTH = false;
    static constexpr bool HAS_START = true;
    static constexpr MoveColsIdx START = MoveColsIdx::START;
};

/* Useful if a specialized move column just extends the existing MoveCols enum, i.e.
   enum class MoreCols {
     LENGTH,
     POINTER,
     OFFSET,
     NEW_VAL1,
     NEW_VAL2,
     NUM_COLS
   };
   
   Just resets the NUM_COLS to the new value */
template <typename Columns, typename BaseColumns>
struct ExtendTraits : public MoveColsTraits<BaseColumns> {
    static constexpr size_t NUM_COLS = static_cast<size_t>(Columns::NUM_COLS);
};

/* Allows for switching between Length and Start columns, if passed one type but want relative or absolute positions.
   Specialized move columns should also implement this
*/
template<typename BaseColumns>
struct ColumnSwitcher;

template<>
struct ColumnSwitcher<MoveCols> {
    using Relative = MoveCols;
    using Absolute = MoveColsIdx;
};

template<>
struct ColumnSwitcher<MoveColsIdx> {
    using Relative = MoveCols;
    using Absolute = MoveColsIdx;
};

// Helper alias for easy switching
template<typename BaseColumns, bool UseAbsolutePositions>
using SwitchColumns = std::conditional_t<UseAbsolutePositions,
    typename ColumnSwitcher<BaseColumns>::Absolute,
    typename ColumnSwitcher<BaseColumns>::Relative>;

/* This was done using generative code... Advanced and not sure how to do it better */
// ADL detector: is C an extended enum (has friend runperm_parent(C))?
template<class C, class = void>
struct is_runperm_extended : std::false_type {};

template<class C>
struct is_runperm_extended<C, std::void_t<decltype(runperm_parent(std::declval<C>()))>>
    : std::true_type {};

// Get the parent RunDataColumns<RD,Base> for extended enums via ADL
template<class C>
using runperm_parent_t = std::remove_pointer_t<decltype(runperm_parent(std::declval<C>()))>;

// Default resolver: plain traits
template<class C, bool = is_runperm_extended<C>::value>
struct ResolveColsTraits { using type = MoveColsTraits<C>; };

// Extended enum: inherit base traits, remap enum constants to C, update NUM_COLS
template<class C>
struct ResolveColsTraits<C, true> {
    using Parent   = runperm_parent_t<C>;           // RunDataColumns<RD, Base>
    using BaseCols = typename Parent::Base;
    using BaseT    = MoveColsTraits<BaseCols>;

    struct type : BaseT {
        using Columns = C;

        static constexpr size_t NUM_COLS = static_cast<size_t>(Columns::NUM_COLS);

        static constexpr Columns PRIMARY = static_cast<Columns>(BaseT::PRIMARY);
        static constexpr Columns POINTER = static_cast<Columns>(BaseT::POINTER);
        static constexpr Columns OFFSET  = static_cast<Columns>(BaseT::OFFSET);

        static constexpr bool HAS_LENGTH = BaseT::HAS_LENGTH;
        static constexpr bool HAS_START  = BaseT::HAS_START;

        using Position = typename BaseT::Position;
        static constexpr bool RELATIVE = BaseT::RELATIVE;
        static constexpr bool ABSOLUTE = BaseT::ABSOLUTE;
    };
};

#endif /* end of include guard: _MOVE_COLUMNS_HPP */