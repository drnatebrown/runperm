// Unit tests for column type traits used by move_structure, permutation, and RLBWT.
// These are simple assert-based tests, no external framework.

#include "orbit/internal/move/move_columns.hpp"
#include "orbit/internal/perm/data_columns.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_columns.hpp"

#include <cassert>
#include <iostream>

using std::size_t;

using namespace orbit;
using namespace orbit::rlbwt;

void test_move_columns_basic_traits() {
    // move_columns (relative)
    using MoveRelTraits = move_cols_traits<move_columns>;
    static_assert(MoveRelTraits::RELATIVE, "move_columns should use relative positions");
    static_assert(MoveRelTraits::PRIMARY == move_columns::LENGTH, "PRIMARY must alias LENGTH");
    static_assert(MoveRelTraits::LENGTH == move_columns::LENGTH, "LENGTH alias mismatch");
    static_assert(MoveRelTraits::POINTER == move_columns::POINTER, "POINTER alias mismatch");
    static_assert(MoveRelTraits::OFFSET == move_columns::OFFSET, "OFFSET alias mismatch");
    static_assert(
        MoveRelTraits::NUM_COLS == static_cast<size_t>(move_columns::COUNT),
        "NUM_COLS must match move_columns::COUNT"
    );

    MoveRelTraits::position rel_pos{};
    rel_pos.interval = 1;
    rel_pos.offset = 2;
    MoveRelTraits::position rel_pos2{};
    rel_pos2.interval = 1;
    rel_pos2.offset = 2;
    MoveRelTraits::position rel_pos3{};
    rel_pos3.interval = 2;
    rel_pos3.offset = 0;

    assert(rel_pos == rel_pos2);
    assert(rel_pos != rel_pos3);

    // move_columns_idx (absolute)
    using MoveAbsTraits = move_cols_traits<move_columns_idx>;
    static_assert(!MoveAbsTraits::RELATIVE, "move_columns_idx should use absolute positions");
    static_assert(MoveAbsTraits::PRIMARY == move_columns_idx::START, "PRIMARY must alias START");
    static_assert(MoveAbsTraits::START == move_columns_idx::START, "START alias mismatch");
    static_assert(MoveAbsTraits::POINTER == move_columns_idx::POINTER, "POINTER alias mismatch");
    static_assert(MoveAbsTraits::OFFSET == move_columns_idx::OFFSET, "OFFSET alias mismatch");
    static_assert(
        MoveAbsTraits::NUM_COLS == static_cast<size_t>(move_columns_idx::COUNT),
        "NUM_COLS must match move_columns_idx::COUNT"
    );

    MoveAbsTraits::position abs_pos{};
    abs_pos.interval = 3;
    abs_pos.offset = 4;
    abs_pos.idx = 7;
    MoveAbsTraits::position abs_pos2{};
    abs_pos2.interval = 3;
    abs_pos2.offset = 4;
    abs_pos2.idx = 10; // idx is not compared in operator==
    MoveAbsTraits::position abs_pos3{};
    abs_pos3.interval = 0;
    abs_pos3.offset = 0;
    abs_pos3.idx = 0;

    assert(abs_pos == abs_pos2);
    assert(abs_pos != abs_pos3);
}

void test_column_switcher_and_switch_columns() {
    // ColumnSwitcher for move_columns
    using switch_from_rel = column_switcher<move_columns>;
    static_assert(std::is_same<switch_from_rel::relative, move_columns>::value,
                  "column_switcher<move_columns>::relative must be move_columns");
    static_assert(std::is_same<switch_from_rel::absolute, move_columns_idx>::value,
                  "column_switcher<move_columns>::absolute must be move_columns_idx");

    // ColumnSwitcher for move_columns_idx
    using switch_from_abs = column_switcher<move_columns_idx>;
    static_assert(std::is_same<switch_from_abs::relative, move_columns>::value,
                  "column_switcher<move_columns_idx>::relative must be move_columns");
    static_assert(std::is_same<switch_from_abs::absolute, move_columns_idx>::value,
                  "column_switcher<move_columns_idx>::absolute must be move_columns_idx");

    // SwitchColumns alias
    using from_rel_abs = switch_columns<move_columns, true>;
    using from_rel_rel = switch_columns<move_columns, false>;
    using from_abs_abs = switch_columns<move_columns_idx, true>;
    using from_abs_rel = switch_columns<move_columns_idx, false>;

    static_assert(std::is_same<from_rel_abs, move_columns_idx>::value,
                  "switch_columns<move_columns,true> must yield move_columns_idx");
    static_assert(std::is_same<from_rel_rel, move_columns>::value,
                  "switch_columns<move_columns,false> must yield move_columns");
    static_assert(std::is_same<from_abs_abs, move_columns_idx>::value,
                  "switch_columns<move_columns_idx,true> must yield move_columns_idx");
    static_assert(std::is_same<from_abs_rel, move_columns>::value,
                  "switch_columns<move_columns_idx,false> must yield move_columns");
}

void test_rlbwt_columns_traits_and_switcher() {
    // Relative RLBWT columns
    using rlbwt_rel_traits = move_cols_traits<rlbwt_columns>;
    static_assert(rlbwt_rel_traits::RELATIVE, "rlbwt_columns should be relative");
    static_assert(rlbwt_rel_traits::PRIMARY == rlbwt_columns::LENGTH, "primary must alias length");
    static_assert(rlbwt_rel_traits::LENGTH == rlbwt_columns::LENGTH, "length alias mismatch");
    static_assert(rlbwt_rel_traits::POINTER == rlbwt_columns::POINTER, "pointer alias mismatch");
    static_assert(rlbwt_rel_traits::OFFSET == rlbwt_columns::OFFSET, "offset alias mismatch");
    static_assert(rlbwt_rel_traits::CHARACTER == rlbwt_columns::CHARACTER, "character alias mismatch");
    static_assert(
        rlbwt_rel_traits::NUM_COLS == static_cast<size_t>(rlbwt_columns::COUNT),
        "num_cols must match rlbwt_columns::COUNT"
    );

    rlbwt_rel_traits::position rel_pos{};
    rel_pos.interval = 5;
    rel_pos.offset = 6;
    rlbwt_rel_traits::position rel_pos2{};
    rel_pos2.interval = 5;
    rel_pos2.offset = 6;
    rlbwt_rel_traits::position rel_pos3{};
    rel_pos3.interval = 5;
    rel_pos3.offset = 0;

    assert(rel_pos == rel_pos2);
    assert(rel_pos != rel_pos3);

    // Absolute RLBWT columns
    using rlbwt_abs_traits = move_cols_traits<rlbwt_columns_idx>;
    static_assert(!rlbwt_abs_traits::RELATIVE, "rlbwt_columns_idx should be absolute");
    static_assert(rlbwt_abs_traits::PRIMARY == rlbwt_columns_idx::START, "primary must alias start");
    static_assert(rlbwt_abs_traits::START == rlbwt_columns_idx::START, "start alias mismatch");
    static_assert(rlbwt_abs_traits::POINTER == rlbwt_columns_idx::POINTER, "pointer alias mismatch");
    static_assert(rlbwt_abs_traits::OFFSET == rlbwt_columns_idx::OFFSET, "offset alias mismatch");
    static_assert(rlbwt_abs_traits::CHARACTER == rlbwt_columns_idx::CHARACTER, "character alias mismatch");
    static_assert(
        rlbwt_abs_traits::NUM_COLS == static_cast<size_t>(rlbwt_columns_idx::COUNT),
        "num_cols must match rlbwt_columns_idx::COUNT"
    );

    rlbwt_abs_traits::position abs_pos{};
    abs_pos.interval = 1;
    abs_pos.offset = 0;
    abs_pos.idx = 10;
    rlbwt_abs_traits::position abs_pos2{};
    abs_pos2.interval = 1;
    abs_pos2.offset = 0;
    abs_pos2.idx = 11;
    rlbwt_abs_traits::position abs_pos3{};
    abs_pos3.interval = 2;
    abs_pos3.offset = 0;
    abs_pos3.idx = 0;

    assert(abs_pos == abs_pos2);
    assert(abs_pos != abs_pos3);

    // ColumnSwitcher specializations for RLBWT columns
    using switch_rlbwt_rel = column_switcher<rlbwt_columns>;
    using switch_rlbwt_abs = column_switcher<rlbwt_columns_idx>;

    static_assert(std::is_same<switch_rlbwt_rel::relative, rlbwt_columns>::value,
                  "column_switcher<rlbwt_columns>::relative must be rlbwt_columns");
    static_assert(std::is_same<switch_rlbwt_rel::absolute, rlbwt_columns_idx>::value,
                  "column_switcher<rlbwt_columns_idx>::absolute must be rlbwt_columns_idx");
    static_assert(std::is_same<switch_rlbwt_abs::relative, rlbwt_columns>::value,
                  "column_switcher<rlbwt_columns_idx>::relative must be rlbwt_columns");
    static_assert(std::is_same<switch_rlbwt_abs::absolute, rlbwt_columns_idx>::value,
                  "column_switcher<rlbwt_columns_idx>::absolute must be rlbwt_columns_idx");
}

// Simple RunColsWrapper usage: ensure extended columns are resolved and indexed correctly.
void test_run_columns_wrapper_and_resolve_traits() {
    // Define a tiny run-data enum.
    enum class my_run_cols {
        A,
        B,
        COUNT
    };

    using wrapper = data_columns_wrapper<my_run_cols, move_columns>;
    static_assert(wrapper::NUM_BASE_COLS == move_cols_traits<move_columns>::NUM_COLS,
                  "num_base_cols must match move_cols_traits<move_columns>::num_cols");
    static_assert(wrapper::NUM_FIELDS == static_cast<size_t>(my_run_cols::COUNT),
                  "num_fields must match my_run_cols::COUNT");
    static_assert(wrapper::NUM_COLS == wrapper::NUM_BASE_COLS + wrapper::NUM_FIELDS,
                  "num_cols must equal base + fields");

    // The nested enum E can be used as a columns type with ResolveColsTraits.
    using columns_e = wrapper::e;
    using traits_e = resolve_cols_traits<columns_e, true>::type;

    static_assert(traits_e::RELATIVE == move_cols_traits<move_columns>::RELATIVE,
                  "Extended traits must inherit RELATIVE from base");
    static_assert(traits_e::NUM_COLS == wrapper::NUM_COLS,
                  "Extended traits num_cols must match wrapper num_cols");

    // Check run_column mapping from run-data column index into the extended enum.
    columns_e a_col = traits_e::template data_column<my_run_cols::A>();
    columns_e b_col = traits_e::template data_column<my_run_cols::B>();

    size_t a_idx = static_cast<size_t>(a_col);
    size_t b_idx = static_cast<size_t>(b_col);

    assert(a_idx == wrapper::NUM_BASE_COLS);
    assert(b_idx == wrapper::NUM_BASE_COLS + 1);
}

int main() {
    test_move_columns_basic_traits();
    test_column_switcher_and_switch_columns();
    test_rlbwt_columns_traits_and_switcher();
    test_run_columns_wrapper_and_resolve_traits();

    std::cout << "columns tests passed" << std::endl;
    return 0;
}
