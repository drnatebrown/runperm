// Integration-style tests for invertible move structure navigation semantics.

#include "orbit/internal/move/move_structure_impl.hpp"
#include "orbit/interval_encoding.hpp"
#include "orbit/move_structure.hpp"

#include <cassert>
#include <iostream>
#include <vector>

using std::size_t;
using std::vector;

using namespace orbit;

static const vector<ulint> kRunnyPerm = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};

template <typename MS>
static ulint global_index_relative(const MS& ms, typename MS::position pos) {
    ulint idx = 0;
    for (size_t i = 0; i < pos.interval; ++i) {
        idx += ms.get_length(i);
    }
    return idx + pos.offset;
}

template <typename MS>
static typename MS::position position_from_index_relative(const MS& ms, ulint idx) {
    ulint prefix = 0;
    for (ulint interval = 0; interval < ms.intervals(); ++interval) {
        const ulint len = ms.get_length(interval);
        if (idx < prefix + len) {
            return {interval, idx - prefix};
        }
        prefix += len;
    }
    assert(false && "index out of range");
    return {};
}

static void test_invertible_splitting_preserves_fwd_inv_semantics() {
    const auto enc_base = invertible_interval_encoding::from_permutation(kRunnyPerm, NO_SPLITTING);
    const auto enc_capped = invertible_interval_encoding::from_permutation(kRunnyPerm, ONLY_LENGTH_CAPPING);

    invertible_structure_vec ms_base(enc_base);
    invertible_structure_vec ms_capped(enc_capped);
    const vector<ulint> inverse = get_inverse_permutation(kRunnyPerm);

    assert(ms_base.domain() == kRunnyPerm.size());
    assert(ms_capped.domain() == kRunnyPerm.size());

    for (ulint idx = 0; idx < kRunnyPerm.size(); ++idx) {
        auto pos_base = position_from_index_relative(ms_base, idx);
        auto pos_cap = position_from_index_relative(ms_capped, idx);

        pos_base = ms_base.move_fwd(pos_base);
        pos_cap = ms_capped.move_fwd(pos_cap);
        assert(global_index_relative(ms_base, pos_base) == global_index_relative(ms_capped, pos_cap));
        assert(global_index_relative(ms_base, pos_base) == kRunnyPerm[idx]);

        pos_base = position_from_index_relative(ms_base, idx);
        pos_cap = position_from_index_relative(ms_capped, idx);
        pos_base = ms_base.move_inv(pos_base);
        pos_cap = ms_capped.move_inv(pos_cap);
        assert(global_index_relative(ms_base, pos_base) == global_index_relative(ms_capped, pos_cap));
        assert(global_index_relative(ms_base, pos_base) == inverse[idx]);
    }
}

static void test_invertible_absolute_splitting_preserves_semantics() {
    const auto enc_base = invertible_interval_encoding::from_permutation(kRunnyPerm, NO_SPLITTING);
    const auto enc_capped = invertible_interval_encoding::from_permutation(kRunnyPerm, ONLY_LENGTH_CAPPING);

    invertible_structure_vec_idx ms_base(enc_base);
    invertible_structure_vec_idx ms_capped(enc_capped);

    for (ulint idx = 0; idx < kRunnyPerm.size(); ++idx) {
        typename invertible_structure_vec_idx::position pos_base{};
        pos_base.idx = idx;
        ulint prefix = 0;
        for (ulint interval = 0; interval < ms_base.intervals(); ++interval) {
            const ulint len = ms_base.get_length(interval);
            if (idx < prefix + len) {
                pos_base.interval = interval;
                pos_base.offset = idx - prefix;
                break;
            }
            prefix += len;
        }

        typename invertible_structure_vec_idx::position pos_cap = pos_base;
        prefix = 0;
        for (ulint interval = 0; interval < ms_capped.intervals(); ++interval) {
            const ulint len = ms_capped.get_length(interval);
            if (idx < prefix + len) {
                pos_cap.interval = interval;
                pos_cap.offset = idx - prefix;
                pos_cap.idx = ms_capped.get_start(interval) + pos_cap.offset;
                break;
            }
            prefix += len;
        }

        pos_base = ms_base.move_fwd(pos_base);
        pos_cap = ms_capped.move_fwd(pos_cap);
        assert(pos_base.idx == pos_cap.idx);
        assert(pos_base.idx == kRunnyPerm[idx]);
    }
}

int main() {
    test_invertible_splitting_preserves_fwd_inv_semantics();
    test_invertible_absolute_splitting_preserves_semantics();

    std::cout << "invertible_structure integration tests passed" << std::endl;
    return 0;
}
