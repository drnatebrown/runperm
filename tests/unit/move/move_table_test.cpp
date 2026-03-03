// Unit tests for MoveRow and MoveTable.
// These are simple assert-based tests, no external framework.

#include "internal/move/move_row.hpp"
#include "internal/move/move_table.hpp"

#include <array>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using std::array;
using std::size_t;
using std::vector;

void test_move_row_relative_basic() {
    using Row = MoveRow<MoveCols>;

    Row row;

    ulint primary = static_cast<ulint>(123) & MASK(Row::RowTraits::PRIMARY_BITS);
    ulint pointer = static_cast<ulint>(456) & MASK(Row::RowTraits::POINTER_BITS);
    ulint offset  = static_cast<ulint>(789) & MASK(Row::RowTraits::OFFSET_BITS);

    row.set<MoveCols::LENGTH>(primary);
    row.set<MoveCols::POINTER>(pointer);
    row.set<MoveCols::OFFSET>(offset);

    assert(row.get<MoveCols::LENGTH>() == primary);
    assert(row.get<MoveCols::POINTER>() == pointer);
    assert(row.get<MoveCols::OFFSET>() == offset);

    constexpr size_t NumCols = static_cast<size_t>(MoveCols::COUNT);
    array<ulint, NumCols> values{};
    values[static_cast<size_t>(MoveCols::LENGTH)]  = primary;
    values[static_cast<size_t>(MoveCols::POINTER)] = pointer;
    values[static_cast<size_t>(MoveCols::OFFSET)]  = offset;

    Row row_from_array(values);
    auto roundtrip = row_from_array.get();
    assert(roundtrip[static_cast<size_t>(MoveCols::LENGTH)]  == primary);
    assert(roundtrip[static_cast<size_t>(MoveCols::POINTER)] == pointer);
    assert(roundtrip[static_cast<size_t>(MoveCols::OFFSET)]  == offset);
}

void test_move_row_absolute_basic() {
    using Row = MoveRow<MoveColsIdx>;

    Row row;

    ulint start   = static_cast<ulint>(42) & MASK(Row::RowTraits::PRIMARY_BITS);
    ulint pointer = static_cast<ulint>(17) & MASK(Row::RowTraits::POINTER_BITS);
    ulint offset  = static_cast<ulint>(9)  & MASK(Row::RowTraits::OFFSET_BITS);

    row.set<MoveColsIdx::START>(start);
    row.set<MoveColsIdx::POINTER>(pointer);
    row.set<MoveColsIdx::OFFSET>(offset);

    assert(row.get<MoveColsIdx::START>() == start);
    assert(row.get<MoveColsIdx::POINTER>() == pointer);
    assert(row.get<MoveColsIdx::OFFSET>() == offset);

    constexpr size_t NumCols = static_cast<size_t>(MoveColsIdx::COUNT);
    array<ulint, NumCols> values{};
    values[static_cast<size_t>(MoveColsIdx::START)]   = start;
    values[static_cast<size_t>(MoveColsIdx::POINTER)] = pointer;
    values[static_cast<size_t>(MoveColsIdx::OFFSET)]  = offset;

    Row row_from_array(values);
    auto roundtrip = row_from_array.get();
    assert(roundtrip[static_cast<size_t>(MoveColsIdx::START)]   == start);
    assert(roundtrip[static_cast<size_t>(MoveColsIdx::POINTER)] == pointer);
    assert(roundtrip[static_cast<size_t>(MoveColsIdx::OFFSET)]  == offset);
}

void test_move_row_assert_widths() {
    {
        using Row = MoveRow<MoveCols>;
        constexpr size_t NumCols = static_cast<size_t>(MoveCols::COUNT);
        array<uchar, NumCols> widths{};

        widths[static_cast<size_t>(MoveCols::LENGTH)] =
            static_cast<uchar>(Row::RowTraits::PRIMARY_BITS);
        widths[static_cast<size_t>(MoveCols::POINTER)] =
            static_cast<uchar>(Row::RowTraits::POINTER_BITS);
        widths[static_cast<size_t>(MoveCols::OFFSET)] =
            static_cast<uchar>(Row::RowTraits::OFFSET_BITS);

        Row::assert_widths(widths);
    }

    {
        using Row = MoveRow<MoveColsIdx>;
        constexpr size_t NumCols = static_cast<size_t>(MoveColsIdx::COUNT);
        array<uchar, NumCols> widths{};

        widths[static_cast<size_t>(MoveColsIdx::START)] =
            static_cast<uchar>(Row::RowTraits::PRIMARY_BITS);
        widths[static_cast<size_t>(MoveColsIdx::POINTER)] =
            static_cast<uchar>(Row::RowTraits::POINTER_BITS);
        widths[static_cast<size_t>(MoveColsIdx::OFFSET)] =
            static_cast<uchar>(Row::RowTraits::OFFSET_BITS);

        Row::assert_widths(widths);
    }
}

void test_move_table_from_packed_vector_relative() {
    using Columns = MoveCols;
    using Table = MoveTable<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();

    constexpr size_t NumCols = Table::NumCols;
    const size_t rows = 32;

    PackedVector<Columns> vec(rows, widths);

    for (size_t i = 0; i < rows; ++i) {
        array<ulint, NumCols> values{};
        values[static_cast<size_t>(Columns::LENGTH)] =
            static_cast<ulint>((i + 1) & MASK(widths[static_cast<size_t>(Columns::LENGTH)]));
        values[static_cast<size_t>(Columns::POINTER)] =
            static_cast<ulint>((i * 7) & MASK(widths[static_cast<size_t>(Columns::POINTER)]));
        values[static_cast<size_t>(Columns::OFFSET)] =
            static_cast<ulint>((i * 11) & MASK(widths[static_cast<size_t>(Columns::OFFSET)]));
        vec.set_row(i, values);
    }

    Table table(std::move(vec));
    assert(table.size() == rows);

    for (size_t i = 0; i < rows; ++i) {
        auto row = table.get_row(i);

        ulint expected_length =
            static_cast<ulint>((i + 1) & MASK(widths[static_cast<size_t>(Columns::LENGTH)]));
        ulint expected_pointer =
            static_cast<ulint>((i * 7) & MASK(widths[static_cast<size_t>(Columns::POINTER)]));
        ulint expected_offset =
            static_cast<ulint>((i * 11) & MASK(widths[static_cast<size_t>(Columns::OFFSET)]));

        assert(row[static_cast<size_t>(Columns::LENGTH)] == expected_length);
        assert(row[static_cast<size_t>(Columns::POINTER)] == expected_pointer);
        assert(row[static_cast<size_t>(Columns::OFFSET)] == expected_offset);

        assert(table.get_length(i) == expected_length);
        assert(table.get_pointer(i) == expected_pointer);
        assert(table.get_offset(i) == expected_offset);
        assert(table.get_primary(i) == expected_length);
    }
}

void test_move_table_interface_relative_and_absolute() {
    {
        using Table = MoveTable<MoveCols>;
        using Row = Table::Row;

        Table table;
        table.table = std::vector<Row>(1);

        // For relative columns, primary corresponds to interval length.
        ulint start = 100;
        ulint length = 5;
        table.set_primary(0, start, length);
        table.set_pointer(0, 7);
        table.set_offset(0, 3);

        assert(table.get_primary(0) == length);
        assert(table.get_length(0) == length);
        assert(table.get_pointer(0) == 7);
        assert(table.get_offset(0) == 3);
    }

    {
        using Table = MoveTable<MoveColsIdx>;
        using Row = Table::Row;

        Table table;
        table.table = std::vector<Row>(1);

        // For absolute columns, primary corresponds to the start index.
        ulint start = 42;
        ulint length = 10;
        table.set_primary(0, start, length);
        assert(table.get_primary(0) == start);
        assert(table.get_start(0) == start);

        table.set_start(0, 100);
        assert(table.get_start(0) == 100);
    }
}

void test_move_table_serialize_roundtrip() {
    using Columns = MoveCols;
    using Table = MoveTable<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();
    constexpr size_t NumCols = Table::NumCols;

    const size_t rows = 16;
    PackedVector<Columns> vec(rows, widths);

    for (size_t i = 0; i < rows; ++i) {
        array<ulint, NumCols> values{};
        values[static_cast<size_t>(Columns::LENGTH)] = static_cast<ulint>(i + 1);
        values[static_cast<size_t>(Columns::POINTER)] = static_cast<ulint>(i * 3);
        values[static_cast<size_t>(Columns::OFFSET)] = static_cast<ulint>(i * 2);
        vec.set_row(i, values);
    }

    Table table(std::move(vec));

    std::stringstream ss;
    size_t bytes_written = table.serialize(ss);
    assert(bytes_written > 0);

    Table loaded;
    loaded.load(ss);

    assert(loaded.size() == table.size());
    for (size_t i = 0; i < rows; ++i) {
        auto original = table.get_row(i);
        auto roundtrip = loaded.get_row(i);
        assert(original == roundtrip);
    }
}

void test_move_table_bits_needed() {
    using Columns = MoveCols;
    using Table = MoveTable<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();

    const size_t rows = 10;
    size_t expected_bits = BYTES_TO_BITS(sizeof(typename Table::Row)) * rows;
    size_t computed_bits = Table::bits_needed(rows, widths);
    assert(computed_bits == expected_bits);
}

int main() {
    test_move_row_relative_basic();
    test_move_row_absolute_basic();
    test_move_row_assert_widths();

    test_move_table_from_packed_vector_relative();
    test_move_table_interface_relative_and_absolute();
    test_move_table_serialize_roundtrip();
    test_move_table_bits_needed();

    std::cout << "move_row/move_table tests passed" << std::endl;
    return 0;
}

