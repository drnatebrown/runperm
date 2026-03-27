// Unit tests for common utilities/macros in `common.hpp`.
// These are simple assert-based tests, consistent with the rest of the suite.

#include "orbit/common.hpp"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

void test_define_run_cols_macro() {
    // Exercise DEFINE_ENUM_CLASS_WITH_COUNT and enum helpers from common.hpp.
    DEFINE_ORBIT_COLUMNS(run_cols, A, B);

    // COUNT is appended as the last enumerator, but the total size is the number of fields provided.
    static_assert(num_columns<run_cols>() == 2, "num_columns should match number of fields provided");

    columns_tuple<run_cols> tuple{};
    tuple[to_index(run_cols::A)] = 7;
    tuple[to_index(run_cols::B)] = 11;

    assert(tuple[to_index(run_cols::A)] == 7);
    assert(tuple[to_index(run_cols::B)] == 11);
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

enum class test_enum {
    A,
    B,
    C,
    COUNT
};

void test_enum_helpers() {
    // to_index should map enum values to their underlying indices.
    assert(to_index(test_enum::A) == 0);
    assert(to_index(test_enum::B) == 1);
    assert(to_index(test_enum::C) == 2);

    // enum_count uses the COUNT sentinel.
    assert(num_columns<test_enum>() == 3);

    // DataTuple builds a fixed-size array consistent with enum_count.
    columns_tuple<test_enum> tuple{};
    static_assert(std::tuple_size<decltype(tuple)>::value == 3,
                  "columns_tuple<test_enum> must have size 3");
    tuple[to_index(test_enum::A)] = 10;
    tuple[to_index(test_enum::B)] = 20;
    tuple[to_index(test_enum::C)] = 30;

    ulint sum = 0;
    for (ulint v : tuple) {
        sum += v;
    }
    assert(sum == 60);
}

void test_macros_sanity() {
    // bytes_to_bits / bits_to_bytes round-trip for multiples of a byte.
    static_assert(bytes_to_bits(1) == 8, "1 byte must be 8 bits");
    static_assert(bytes_to_bits(4) == 32, "4 bytes must be 32 bits");
    static_assert(bits_to_bytes(8) == 1, "8 bits must be 1 byte");
    static_assert(bits_to_bytes(32) == 4, "32 bits must be 4 bytes");

    // num_bits_type agrees with sizeof.
    static_assert(num_bits_type(ulint) == bytes_to_bits(sizeof(ulint)),
                  "num_bits_type must be sizeof(type) * 8");

    // ceil_div basic behaviour.
    assert(ceil_div(0, 1) == 0);
    assert(ceil_div(1, 1) == 1);
    assert(ceil_div(5, 2) == 3);
    assert(ceil_div(6, 2) == 3);
    assert(ceil_div(7, 2) == 4);

    // pow2 / max_val / mask basic checks.
    static_assert(pow2(0) == 1ULL, "2^0 must be 1");
    static_assert(pow2(3) == 8ULL, "2^3 must be 8");

    static_assert(max_val(1) == 1ULL, "max_val(1) = 1");
    static_assert(max_val(3) == 7ULL, "max_val(3) = 7");

    static_assert(mask(4) == 0xFULL, "mask(4) must be 0b1111");
}

int main() {
    test_bit_width_basic();
    test_define_run_cols_macro();
    test_enum_helpers();
    test_macros_sanity();

    std::cout << "orbit/common tests passed" << std::endl;
    return 0;
}
