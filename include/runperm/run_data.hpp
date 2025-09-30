#ifndef _RUN_DATA_HPP
#define _RUN_DATA_HPP

#include "move/move_columns.hpp"

template<typename RunData, typename BaseColumns = MoveCols>
struct RunDataColumns {
    using RD = RunData;
    using Base = BaseColumns;
    using BaseTraits = MoveColsTraits<BaseColumns>;
    // Inherits from some BaseColumns type used by the move structure
    enum class E : size_t {
        // mirror base columns
        PRIMARY  = static_cast<size_t>(BaseTraits::PRIMARY),
        POINTER = static_cast<size_t>(BaseTraits::POINTER),
        OFFSET  = static_cast<size_t>(BaseTraits::OFFSET),
        /* other unlisted base columns here*/
        // total columns = base + user fields
        NUM_COLS = static_cast<size_t>(BaseTraits::NUM_COLS) + static_cast<size_t>(RunData::NUM_COLS)
    };
    
    template<RunData RunDataField>
    static constexpr E field() { 
        return static_cast<E>(static_cast<size_t>(BaseTraits::NUM_COLS) + static_cast<size_t>(RunDataField)); 
    }
    
    static constexpr size_t NUM_BASE_COLS = static_cast<size_t>(BaseTraits::NUM_COLS);
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunData::NUM_COLS);
    static constexpr size_t NUM_COLS = NUM_BASE_COLS + NUM_FIELDS;

    // ADL back-link
    friend constexpr RunDataColumns* runperm_parent(E) { return nullptr; }
};

#endif /* end of include guard: _RUN_DATA_HPP */