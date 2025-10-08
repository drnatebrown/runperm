#ifndef _RUN_DATA_HPP
#define _RUN_DATA_HPP

#include "move/move_columns.hpp"

// Special empty run columns type for when no run data is needed
enum class EmptyRunCols {
    COUNT // Helper to get the number of columns
};

template<typename RunColsType, typename BaseColumns = MoveCols>
struct RunColsWrapper {
    using RunCols = RunColsType;
    using Base = BaseColumns;
    using BaseTraits = MoveColsTraits<BaseColumns>;

    static constexpr size_t NUM_BASE_COLS = static_cast<size_t>(BaseTraits::NUM_COLS);
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunCols::COUNT);
    static constexpr size_t NUM_COLS = NUM_BASE_COLS + NUM_FIELDS;
    
    // "trick" to pass a columns type that has everything needed for the move structure
    // even if it doesn't explicitly have them. Importantly, the count is correct.
    enum class E : size_t {
        /* all columns unlisted here */
        // total columns = base + user fields
        COUNT = NUM_COLS
    };

    // Back link to parent structure, allows E to be used as a columns type
    // with its own traits, inherit from BaseColumns
    friend constexpr RunColsWrapper* runperm_parent(E) { return nullptr; }
};

// Extended enum: inherit base traits, update NUM_COLS, and allow special function to resolve run columns
template<class C>
struct ResolveColsTraits<C, true> {
    using Parent   = runperm_parent_t<C>;           // RunDataColumns<RunData, Base>
    using RunCols = typename Parent::RunCols;
    using BaseCols = typename Parent::Base;
    using BaseTraits    = MoveColsTraits<BaseCols>;

    struct type : BaseTraits {
        using Columns = C;
        static constexpr size_t NUM_COLS = static_cast<size_t>(Parent::NUM_COLS);

        template<RunCols RunDataCol>
        static constexpr C run_column() { 
            return static_cast<C>(Parent::NUM_BASE_COLS + static_cast<size_t>(RunDataCol)); 
        }
    };
};

#endif /* end of include guard: _RUN_DATA_HPP */