// Unit tests for column type traits used by move, runperm, and RLBWT.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_columns.hpp"
#include "internal/runperm/run_columns.hpp"
#include "internal/rlbwt/specializations/rlbwt_columns.hpp"

#include <cassert>
#include <iostream>

using std::size_t;

void test_move_columns_basic_traits() {
    // MoveCols (relative)
    using MoveRelTraits = MoveColsTraits<MoveCols>;
    static_assert(MoveRelTraits::RELATIVE, "MoveCols should use relative positions");
    static_assert(MoveRelTraits::PRIMARY == MoveCols::LENGTH, "PRIMARY must alias LENGTH");
    static_assert(MoveRelTraits::LENGTH == MoveCols::LENGTH, "LENGTH alias mismatch");
    static_assert(MoveRelTraits::POINTER == MoveCols::POINTER, "POINTER alias mismatch");
    static_assert(MoveRelTraits::OFFSET == MoveCols::OFFSET, "OFFSET alias mismatch");
    static_assert(
        MoveRelTraits::NUM_COLS == static_cast<size_t>(MoveCols::COUNT),
        "NUM_COLS must match MoveCols::COUNT"
    );

    MoveRelTraits::Position rel_pos{};
    rel_pos.interval = 1;
    rel_pos.offset = 2;
    MoveRelTraits::Position rel_pos2{};
    rel_pos2.interval = 1;
    rel_pos2.offset = 2;
    MoveRelTraits::Position rel_pos3{};
    rel_pos3.interval = 2;
    rel_pos3.offset = 0;

    assert(rel_pos == rel_pos2);
    assert(rel_pos != rel_pos3);

    // MoveColsIdx (absolute)
    using MoveAbsTraits = MoveColsTraits<MoveColsIdx>;
    static_assert(!MoveAbsTraits::RELATIVE, "MoveColsIdx should use absolute positions");
    static_assert(MoveAbsTraits::PRIMARY == MoveColsIdx::START, "PRIMARY must alias START");
    static_assert(MoveAbsTraits::START == MoveColsIdx::START, "START alias mismatch");
    static_assert(MoveAbsTraits::POINTER == MoveColsIdx::POINTER, "POINTER alias mismatch");
    static_assert(MoveAbsTraits::OFFSET == MoveColsIdx::OFFSET, "OFFSET alias mismatch");
    static_assert(
        MoveAbsTraits::NUM_COLS == static_cast<size_t>(MoveColsIdx::COUNT),
        "NUM_COLS must match MoveColsIdx::COUNT"
    );

    MoveAbsTraits::Position abs_pos{};
    abs_pos.interval = 3;
    abs_pos.offset = 4;
    abs_pos.idx = 7;
    MoveAbsTraits::Position abs_pos2{};
    abs_pos2.interval = 3;
    abs_pos2.offset = 4;
    abs_pos2.idx = 10; // idx is not compared in operator==
    MoveAbsTraits::Position abs_pos3{};
    abs_pos3.interval = 0;
    abs_pos3.offset = 0;
    abs_pos3.idx = 0;

    assert(abs_pos == abs_pos2);
    assert(abs_pos != abs_pos3);
}

void test_column_switcher_and_switch_columns() {
    // ColumnSwitcher for MoveCols
    using SwitchFromRel = ColumnSwitcher<MoveCols>;
    static_assert(std::is_same<SwitchFromRel::Relative, MoveCols>::value,
                  "ColumnSwitcher<MoveCols>::Relative must be MoveCols");
    static_assert(std::is_same<SwitchFromRel::Absolute, MoveColsIdx>::value,
                  "ColumnSwitcher<MoveCols>::Absolute must be MoveColsIdx");

    // ColumnSwitcher for MoveColsIdx
    using SwitchFromAbs = ColumnSwitcher<MoveColsIdx>;
    static_assert(std::is_same<SwitchFromAbs::Relative, MoveCols>::value,
                  "ColumnSwitcher<MoveColsIdx>::Relative must be MoveCols");
    static_assert(std::is_same<SwitchFromAbs::Absolute, MoveColsIdx>::value,
                  "ColumnSwitcher<MoveColsIdx>::Absolute must be MoveColsIdx");

    // SwitchColumns alias
    using FromRelAbs = SwitchColumns<MoveCols, true>;
    using FromRelRel = SwitchColumns<MoveCols, false>;
    using FromAbsAbs = SwitchColumns<MoveColsIdx, true>;
    using FromAbsRel = SwitchColumns<MoveColsIdx, false>;

    static_assert(std::is_same<FromRelAbs, MoveColsIdx>::value,
                  "SwitchColumns<MoveCols,true> must yield MoveColsIdx");
    static_assert(std::is_same<FromRelRel, MoveCols>::value,
                  "SwitchColumns<MoveCols,false> must yield MoveCols");
    static_assert(std::is_same<FromAbsAbs, MoveColsIdx>::value,
                  "SwitchColumns<MoveColsIdx,true> must yield MoveColsIdx");
    static_assert(std::is_same<FromAbsRel, MoveCols>::value,
                  "SwitchColumns<MoveColsIdx,false> must yield MoveCols");
}

void test_rlbwt_columns_traits_and_switcher() {
    // Relative RLBWT columns
    using RlbwtRelTraits = MoveColsTraits<RLBWTCols>;
    static_assert(RlbwtRelTraits::RELATIVE, "RLBWTCols should be relative");
    static_assert(RlbwtRelTraits::PRIMARY == RLBWTCols::LENGTH, "PRIMARY must alias LENGTH");
    static_assert(RlbwtRelTraits::LENGTH == RLBWTCols::LENGTH, "LENGTH alias mismatch");
    static_assert(RlbwtRelTraits::POINTER == RLBWTCols::POINTER, "POINTER alias mismatch");
    static_assert(RlbwtRelTraits::OFFSET == RLBWTCols::OFFSET, "OFFSET alias mismatch");
    static_assert(RlbwtRelTraits::CHARACTER == RLBWTCols::CHARACTER, "CHARACTER alias mismatch");
    static_assert(
        RlbwtRelTraits::NUM_COLS == static_cast<size_t>(RLBWTCols::COUNT),
        "NUM_COLS must match RLBWTCols::COUNT"
    );

    RlbwtRelTraits::Position rel_pos{};
    rel_pos.interval = 5;
    rel_pos.offset = 6;
    RlbwtRelTraits::Position rel_pos2{};
    rel_pos2.interval = 5;
    rel_pos2.offset = 6;
    RlbwtRelTraits::Position rel_pos3{};
    rel_pos3.interval = 5;
    rel_pos3.offset = 0;

    assert(rel_pos == rel_pos2);
    assert(rel_pos != rel_pos3);

    // Absolute RLBWT columns
    using RlbwtAbsTraits = MoveColsTraits<RLBWTColsIdx>;
    static_assert(!RlbwtAbsTraits::RELATIVE, "RLBWTColsIdx should be absolute");
    static_assert(RlbwtAbsTraits::PRIMARY == RLBWTColsIdx::START, "PRIMARY must alias START");
    static_assert(RlbwtAbsTraits::START == RLBWTColsIdx::START, "START alias mismatch");
    static_assert(RlbwtAbsTraits::POINTER == RLBWTColsIdx::POINTER, "POINTER alias mismatch");
    static_assert(RlbwtAbsTraits::OFFSET == RLBWTColsIdx::OFFSET, "OFFSET alias mismatch");
    static_assert(RlbwtAbsTraits::CHARACTER == RLBWTColsIdx::CHARACTER, "CHARACTER alias mismatch");
    static_assert(
        RlbwtAbsTraits::NUM_COLS == static_cast<size_t>(RLBWTColsIdx::COUNT),
        "NUM_COLS must match RLBWTColsIdx::COUNT"
    );

    RlbwtAbsTraits::Position abs_pos{};
    abs_pos.interval = 1;
    abs_pos.offset = 0;
    abs_pos.idx = 10;
    RlbwtAbsTraits::Position abs_pos2{};
    abs_pos2.interval = 1;
    abs_pos2.offset = 0;
    abs_pos2.idx = 11;
    RlbwtAbsTraits::Position abs_pos3{};
    abs_pos3.interval = 2;
    abs_pos3.offset = 0;
    abs_pos3.idx = 0;

    assert(abs_pos == abs_pos2);
    assert(abs_pos != abs_pos3);

    // ColumnSwitcher specializations for RLBWT columns
    using SwitchRlbwtRel = ColumnSwitcher<RLBWTCols>;
    using SwitchRlbwtAbs = ColumnSwitcher<RLBWTColsIdx>;

    static_assert(std::is_same<SwitchRlbwtRel::Relative, RLBWTCols>::value,
                  "ColumnSwitcher<RLBWTCols>::Relative must be RLBWTCols");
    static_assert(std::is_same<SwitchRlbwtRel::Absolute, RLBWTColsIdx>::value,
                  "ColumnSwitcher<RLBWTCols>::Absolute must be RLBWTColsIdx");
    static_assert(std::is_same<SwitchRlbwtAbs::Relative, RLBWTCols>::value,
                  "ColumnSwitcher<RLBWTColsIdx>::Relative must be RLBWTCols");
    static_assert(std::is_same<SwitchRlbwtAbs::Absolute, RLBWTColsIdx>::value,
                  "ColumnSwitcher<RLBWTColsIdx>::Absolute must be RLBWTColsIdx");
}

// Simple RunColsWrapper usage: ensure extended columns are resolved and indexed correctly.
void test_run_columns_wrapper_and_resolve_traits() {
    // Define a tiny run-data enum.
    enum class MyRunCols {
        A,
        B,
        COUNT
    };

    using Wrapper = RunColsWrapper<MyRunCols, MoveCols>;
    static_assert(Wrapper::NUM_BASE_COLS == MoveColsTraits<MoveCols>::NUM_COLS,
                  "NUM_BASE_COLS must match MoveColsTraits<MoveCols>::NUM_COLS");
    static_assert(Wrapper::NUM_FIELDS == static_cast<size_t>(MyRunCols::COUNT),
                  "NUM_FIELDS must match MyRunCols::COUNT");
    static_assert(Wrapper::NUM_COLS == Wrapper::NUM_BASE_COLS + Wrapper::NUM_FIELDS,
                  "NUM_COLS must equal base + fields");

    // The nested enum E can be used as a columns type with ResolveColsTraits.
    using ColumnsE = Wrapper::E;
    using TraitsE = ResolveColsTraits<ColumnsE, true>::type;

    static_assert(TraitsE::RELATIVE == MoveColsTraits<MoveCols>::RELATIVE,
                  "Extended traits must inherit RELATIVE from base");
    static_assert(TraitsE::NUM_COLS == Wrapper::NUM_COLS,
                  "Extended traits NUM_COLS must match wrapper NUM_COLS");

    // Check run_column mapping from run-data column index into the extended enum.
    ColumnsE a_col = TraitsE::template run_column<MyRunCols::A>();
    ColumnsE b_col = TraitsE::template run_column<MyRunCols::B>();

    size_t a_idx = static_cast<size_t>(a_col);
    size_t b_idx = static_cast<size_t>(b_col);

    assert(a_idx == Wrapper::NUM_BASE_COLS);
    assert(b_idx == Wrapper::NUM_BASE_COLS + 1);
}

int main() {
    test_move_columns_basic_traits();
    test_column_switcher_and_switch_columns();
    test_rlbwt_columns_traits_and_switcher();
    test_run_columns_wrapper_and_resolve_traits();

    std::cout << "columns tests passed" << std::endl;
    return 0;
}
