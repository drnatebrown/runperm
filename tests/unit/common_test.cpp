// Unit tests for common utilities/macros in `common.hpp`.
// These are simple assert-based tests, consistent with the rest of the suite.

#include "common.hpp"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using std::size_t;
using std::vector;

void test_define_run_cols_macro() {
    // Exercise DEFINE_ENUM_CLASS_WITH_COUNT and enum helpers from common.hpp.
    DEFINE_COLUMNS(RunCols, A, B);

    // COUNT is appended as the last enumerator, but the total size is the number of fields provided.
    static_assert(num_columns<RunCols>() == 2, "num_columns should match number of fields provided");

    ColumnsTuple<RunCols> tuple{};
    tuple[to_index(RunCols::A)] = 7;
    tuple[to_index(RunCols::B)] = 11;

    assert(tuple[to_index(RunCols::A)] == 7);
    assert(tuple[to_index(RunCols::B)] == 11);
}

void test_bit_width_basic() {
    // By contract, bit_width(0) is 1.
    assert(bit_width(0) == 1);

    // Powers of two.
    assert(bit_width(1) == 1);
    assert(bit_width(2) == 2);
    assert(bit_width(4) == 3);
    assert(bit_width(8) == 4);

    // Non powers of two: bit width is floor(log2(v)) + 1.
    assert(bit_width(3) == 2);
    assert(bit_width(5) == 3);
    assert(bit_width(7) == 3);
    assert(bit_width(9) == 4);
}

enum class TestEnum {
    A,
    B,
    C,
    COUNT
};

void test_enum_helpers() {
    // to_index should map enum values to their underlying indices.
    assert(to_index(TestEnum::A) == 0);
    assert(to_index(TestEnum::B) == 1);
    assert(to_index(TestEnum::C) == 2);

    // enum_count uses the COUNT sentinel.
    assert(num_columns<TestEnum>() == 3);

    // DataTuple builds a fixed-size array consistent with enum_count.
    ColumnsTuple<TestEnum> tuple{};
    static_assert(std::tuple_size<decltype(tuple)>::value == 3,
                  "DataTuple<TestEnum> must have size 3");
    tuple[to_index(TestEnum::A)] = 10;
    tuple[to_index(TestEnum::B)] = 20;
    tuple[to_index(TestEnum::C)] = 30;

    ulint sum = 0;
    for (ulint v : tuple) {
        sum += v;
    }
    assert(sum == 60);
}

void test_bwt_to_rlbwt_basic() {
    {
        // Empty input -> empty output.
        vector<uchar> bwt;
        auto [heads, lengths] = bwt_to_rlbwt(bwt);
        assert(heads.empty());
        assert(lengths.empty());
    }
    {
        // Single run.
        vector<uchar> bwt = {'A', 'A', 'A'};
        auto [heads, lengths] = bwt_to_rlbwt(bwt);
        assert(heads.size() == 1);
        assert(lengths.size() == 1);
        assert(heads[0] == static_cast<uchar>('A'));
        assert(lengths[0] == 3);
    }
    {
        // Multiple runs.
        vector<uchar> bwt = {'A', 'A', 'C', 'C', 'C', 'G', 'G'};
        auto [heads, lengths] = bwt_to_rlbwt(bwt);

        assert(heads.size() == 3);
        assert(lengths.size() == 3);

        assert(heads[0] == static_cast<uchar>('A'));
        assert(lengths[0] == 2);

        assert(heads[1] == static_cast<uchar>('C'));
        assert(lengths[1] == 3);

        assert(heads[2] == static_cast<uchar>('G'));
        assert(lengths[2] == 2);

        // Sum of run lengths equals original length.
        ulint total = 0;
        for (ulint v : lengths) {
            total += v;
        }
        assert(total == bwt.size());
    }
}

void test_macros_sanity() {
    // BYTES_TO_BITS / BITS_TO_BYTES round-trip for multiples of a byte.
    static_assert(BYTES_TO_BITS(1) == 8, "1 byte must be 8 bits");
    static_assert(BYTES_TO_BITS(4) == 32, "4 bytes must be 32 bits");
    static_assert(BITS_TO_BYTES(8) == 1, "8 bits must be 1 byte");
    static_assert(BITS_TO_BYTES(32) == 4, "32 bits must be 4 bytes");

    // NUM_BITS agrees with sizeof.
    static_assert(NUM_BITS(ulint) == BYTES_TO_BITS(sizeof(ulint)),
                  "NUM_BITS must be sizeof(type) * 8");

    // CEIL_DIV basic behaviour.
    assert(CEIL_DIV(0, 1) == 0);
    assert(CEIL_DIV(1, 1) == 1);
    assert(CEIL_DIV(5, 2) == 3);
    assert(CEIL_DIV(6, 2) == 3);
    assert(CEIL_DIV(7, 2) == 4);

    // POW2 / MAX_VAL / MASK basic checks.
    static_assert(POW2(0) == 1ULL, "2^0 must be 1");
    static_assert(POW2(3) == 8ULL, "2^3 must be 8");

    static_assert(MAX_VAL(1) == 1ULL, "MAX_VAL(1) = 1");
    static_assert(MAX_VAL(3) == 7ULL, "MAX_VAL(3) = 7");

    static_assert(MASK(4) == 0xFULL, "MASK(4) must be 0b1111");
}

int main() {
    test_bit_width_basic();
    test_define_run_cols_macro();
    test_enum_helpers();
    test_bwt_to_rlbwt_basic();
    test_macros_sanity();

    std::cout << "common tests passed" << std::endl;
    return 0;
}
