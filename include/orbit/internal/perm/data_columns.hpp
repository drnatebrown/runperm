#ifndef _RUN_COLUMNS_HPP
#define _RUN_COLUMNS_HPP

#include "orbit/internal/move/move_columns.hpp"

namespace orbit {

// Special empty run columns type for when no run data is needed
enum class empty_data_columns {
    COUNT // Helper to get the number of columns
};

template<typename data_columns_t, typename base_columns_t = move_columns>
struct data_columns_wrapper {
    using data_columns = data_columns_t;
    using base_columns = base_columns_t;
    using base_traits = move_cols_traits<base_columns>;

    static constexpr size_t NUM_BASE_COLS = static_cast<size_t>(base_traits::NUM_COLS);
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(data_columns::COUNT);
    static constexpr size_t NUM_COLS = NUM_BASE_COLS + NUM_FIELDS;
    
    // "trick" to pass a columns type that has everything needed for the move structure
    // even if it doesn't explicitly have them. Importantly, the count is correct.
    enum class e : size_t {
        /* all columns unlisted here */
        // total columns = base + user fields
        COUNT = NUM_COLS
    };

    // Back link to parent structure, allows E to be used as a columns type
    // with its own traits, inherit from BaseColumns
    friend constexpr data_columns_wrapper* columns_parent(e) { return nullptr; }
};

// Extended enum: inherit base traits, update NUM_COLS, and allow special function to resolve run columns
template<class c>
struct resolve_cols_traits<c, true> {
    using parent   = columns_parent_t<c>;           // data_columns_wrapper<data_columns, base_columns>
    using data_columns = typename parent::data_columns;
    using base_columns = typename parent::base_columns;
    using base_traits    = move_cols_traits<base_columns>;

    struct type : base_traits {
        using columns = c;
        static constexpr size_t NUM_COLS = static_cast<size_t>(parent::NUM_COLS);

        template<data_columns dc>
        static constexpr c data_column() { 
            return static_cast<c>(parent::NUM_BASE_COLS + static_cast<size_t>(dc)); 
        }
    };
};

} // namespace orbit

#endif /* end of include guard: _RUN_DATA_COLUMNS_HPP */