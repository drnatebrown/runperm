// Unit tests for MoveRow and MoveTable.
// These are simple assert-based tests, no external framework.

#include "orbit/internal/move/move_row.hpp"
#include "orbit/internal/move/move_table.hpp"

#include <array>
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

using std::array;
using std::size_t;
using std::vector;

using namespace orbit;

void test_move_row_relative_basic() {
    using Row = move_row<move_columns>;

    Row row;

    ulint primary = static_cast<ulint>(123) & mask(Row::row_traits::PRIMARY_BITS);
    ulint pointer = static_cast<ulint>(456) & mask(Row::row_traits::POINTER_BITS);
    ulint offset  = static_cast<ulint>(789) & mask(Row::row_traits::OFFSET_BITS);

    row.set<move_columns::LENGTH>(primary);
    row.set<move_columns::POINTER>(pointer);
    row.set<move_columns::OFFSET>(offset);

    assert(row.get<move_columns::LENGTH>() == primary);
    assert(row.get<move_columns::POINTER>() == pointer);
    assert(row.get<move_columns::OFFSET>() == offset);

    constexpr size_t NumCols = static_cast<size_t>(move_columns::COUNT);
    array<ulint, NumCols> values{};
    values[static_cast<size_t>(move_columns::LENGTH)]  = primary;
    values[static_cast<size_t>(move_columns::POINTER)] = pointer;
    values[static_cast<size_t>(move_columns::OFFSET)]  = offset;

    Row row_from_array(values);
    auto roundtrip = row_from_array.get();
    assert(roundtrip[static_cast<size_t>(move_columns::LENGTH)]  == primary);
    assert(roundtrip[static_cast<size_t>(move_columns::POINTER)] == pointer);
    assert(roundtrip[static_cast<size_t>(move_columns::OFFSET)]  == offset);
}

void test_move_row_absolute_basic() {
    using Row = move_row<move_columns_idx>;

    Row row;

    ulint start   = static_cast<ulint>(42) & mask(Row::row_traits::PRIMARY_BITS);
    ulint pointer = static_cast<ulint>(17) & mask(Row::row_traits::POINTER_BITS);
    ulint offset  = static_cast<ulint>(9)  & mask(Row::row_traits::OFFSET_BITS);

    row.set<move_columns_idx::START>(start);
    row.set<move_columns_idx::POINTER>(pointer);
    row.set<move_columns_idx::OFFSET>(offset);

    assert(row.get<move_columns_idx::START>() == start);
    assert(row.get<move_columns_idx::POINTER>() == pointer);
    assert(row.get<move_columns_idx::OFFSET>() == offset);

    constexpr size_t NumCols = static_cast<size_t>(move_columns_idx::COUNT);
    array<ulint, NumCols> values{};
    values[static_cast<size_t>(move_columns_idx::START)]   = start;
    values[static_cast<size_t>(move_columns_idx::POINTER)] = pointer;
    values[static_cast<size_t>(move_columns_idx::OFFSET)]  = offset;

    Row row_from_array(values);
    auto roundtrip = row_from_array.get();
    assert(roundtrip[static_cast<size_t>(move_columns_idx::START)]   == start);
    assert(roundtrip[static_cast<size_t>(move_columns_idx::POINTER)] == pointer);
    assert(roundtrip[static_cast<size_t>(move_columns_idx::OFFSET)]  == offset);
}

void test_move_row_assert_widths() {
    {
        using Row = move_row<move_columns>;
        constexpr size_t NumCols = static_cast<size_t>(move_columns::COUNT);
        array<uchar, NumCols> widths{};

        widths[static_cast<size_t>(move_columns::LENGTH)] =
            static_cast<uchar>(Row::row_traits::PRIMARY_BITS);
        widths[static_cast<size_t>(move_columns::POINTER)] =
            static_cast<uchar>(Row::row_traits::POINTER_BITS);
        widths[static_cast<size_t>(move_columns::OFFSET)] =
            static_cast<uchar>(Row::row_traits::OFFSET_BITS);

        Row::assert_widths(widths);
    }

    {
        using Row = move_row<move_columns_idx>;
        constexpr size_t NumCols = static_cast<size_t>(move_columns_idx::COUNT);
        array<uchar, NumCols> widths{};

        widths[static_cast<size_t>(move_columns_idx::START)] =
            static_cast<uchar>(Row::row_traits::PRIMARY_BITS);
        widths[static_cast<size_t>(move_columns_idx::POINTER)] =
            static_cast<uchar>(Row::row_traits::POINTER_BITS);
        widths[static_cast<size_t>(move_columns_idx::OFFSET)] =
            static_cast<uchar>(Row::row_traits::OFFSET_BITS);

        Row::assert_widths(widths);
    }
}

void test_move_table_from_packed_vector_relative() {
    using Columns = move_columns;
    using Table = move_table<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();

    constexpr size_t NumCols = Table::num_cols;
    const size_t rows = 32;

    packed_vector<Columns> vec(rows, widths);

    for (size_t i = 0; i < rows; ++i) {
        array<ulint, NumCols> values{};
        values[static_cast<size_t>(Columns::LENGTH)] =
            static_cast<ulint>((i + 1) & mask(widths[static_cast<size_t>(Columns::LENGTH)]));
        values[static_cast<size_t>(Columns::POINTER)] =
            static_cast<ulint>((i * 7) & mask(widths[static_cast<size_t>(Columns::POINTER)]));
        values[static_cast<size_t>(Columns::OFFSET)] =
            static_cast<ulint>((i * 11) & mask(widths[static_cast<size_t>(Columns::OFFSET)]));
        vec.set_row(i, values);
    }

    Table table(std::move(vec));
    assert(table.size() == rows);

    for (size_t i = 0; i < rows; ++i) {
        auto row = table.get_row(i);

        ulint expected_length =
            static_cast<ulint>((i + 1) & mask(widths[static_cast<size_t>(Columns::LENGTH)]));
        ulint expected_pointer =
            static_cast<ulint>((i * 7) & mask(widths[static_cast<size_t>(Columns::POINTER)]));
        ulint expected_offset =
            static_cast<ulint>((i * 11) & mask(widths[static_cast<size_t>(Columns::OFFSET)]));

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
        using Table = move_table<move_columns>;
        using Row = Table::row;

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
        using Table = move_table<move_columns_idx>;
        using Row = Table::row;

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
    using Columns = move_columns;
    using Table = move_table<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();
    constexpr size_t NumCols = Table::num_cols;

    const size_t rows = 16;
    packed_vector<Columns> vec(rows, widths);

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

void test_invertible_row_relative_basic() {
    using Row = invertible_row<invertible_columns>;

    Row row;
    const ulint length = static_cast<ulint>(7) & mask(Row::row_traits::PRIMARY_BITS);
    const ulint pointer_fwd = static_cast<ulint>(3) & mask(Row::row_traits::POINTER_FWD_BITS);
    const ulint pointer_inv = static_cast<ulint>(5) & mask(Row::row_traits::POINTER_INV_BITS);

    row.set<invertible_columns::LENGTH>(length);
    row.set<invertible_columns::POINTER_FWD>(pointer_fwd);
    row.set<invertible_columns::POINTER_INV>(pointer_inv);
    row.set<invertible_columns::FWD_INTERVAL>(1);
    row.set<invertible_columns::INV_INTERVAL>(0);

    assert(row.get<invertible_columns::LENGTH>() == length);
    assert(row.get<invertible_columns::POINTER_FWD>() == pointer_fwd);
    assert(row.get<invertible_columns::POINTER_INV>() == pointer_inv);
    assert(row.get<invertible_columns::FWD_INTERVAL>() == 1);
    assert(row.get<invertible_columns::INV_INTERVAL>() == 0);
}

void test_move_table_bits_needed() {
    using Columns = move_columns;
    using Table = move_table<Columns>;

    Table tmp;
    const auto &widths = tmp.get_widths();

    const size_t rows = 10;
    size_t expected_bits = bytes_to_bits(sizeof(typename Table::row)) * rows;
    size_t computed_bits = Table::bits_needed(rows, widths);
    assert(computed_bits == expected_bits);
}

int main() {
    test_move_row_relative_basic();
    test_move_row_absolute_basic();
    test_move_row_assert_widths();
    test_invertible_row_relative_basic();

    test_move_table_from_packed_vector_relative();
    test_move_table_interface_relative_and_absolute();
    test_move_table_serialize_roundtrip();
    test_move_table_bits_needed();

    std::cout << "move_row/move_table tests passed" << std::endl;
    return 0;
}
