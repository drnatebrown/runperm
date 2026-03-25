// Tiny unit tests for phi / phi_inv helpers and wrappers.
// Focus: structural sanity and wrapper behaviour, not full SA reconstruction.

#include "orbit/rlbwt.hpp"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using namespace orbit;
using namespace orbit::rlbwt;

using std::size_t;
using std::vector;

void test_phi_phi_inv_structure_on_small_rlbwt() {
    // Use the same small example as in the integration test, which is known to be valid.
    vector<uchar> heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> lens  =       { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain;
    ulint max_length;
    auto [phi_lengths, phi_perm] = rlbwt_to_phi_images(heads, lens, &phi_domain, &max_length);
    size_t inv_domain;
    ulint max_length_inv;
    auto [inv_lengths, inv_perm] = rlbwt_to_phi_inv_images(heads, lens, &inv_domain, &max_length_inv);

    assert(phi_domain == 27);
    assert(inv_domain == 27);
    assert(phi_lengths.size() == phi_perm.size());
    assert(inv_lengths.size() == inv_perm.size());

    auto sum_phi = std::accumulate(phi_lengths.begin(), phi_lengths.end(), static_cast<ulint>(0));
    auto sum_inv = std::accumulate(inv_lengths.begin(), inv_lengths.end(), static_cast<ulint>(0));
    assert(sum_phi == phi_domain);
    assert(sum_inv == inv_domain);

    // All SA samples must be within domain.
    for (ulint v : phi_perm) {
        assert(v < phi_domain);
    }
    for (ulint v : inv_perm) {
        assert(v < inv_domain);
    }
}

void test_runperm_phi_phi_inv_wrapper_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain;
    ulint max_length;
    auto [phi_lengths, phi_perm] = rlbwt_to_phi_images(bwt_heads, bwt_run_lengths, &phi_domain, &max_length);
    size_t inv_domain;
    ulint max_length_inv;
    auto [inv_lengths, inv_perm] = rlbwt_to_phi_inv_images(bwt_heads, bwt_run_lengths, &inv_domain, &max_length_inv);

    assert(phi_domain == inv_domain);

    enum class RunCols {
        V,
        COUNT
    };
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunCols::COUNT);
    vector<std::array<ulint, NUM_FIELDS>> run_data_phi(phi_lengths.size());
    vector<std::array<ulint, NUM_FIELDS>> run_data_inv(inv_lengths.size());
    for (size_t i = 0; i < run_data_phi.size(); ++i) {
        run_data_phi[i][0] = static_cast<ulint>(i);
    }
    for (size_t i = 0; i < run_data_inv.size(); ++i) {
        run_data_inv[i][0] = static_cast<ulint>(i * 2);
    }

    runperm_phi<RunCols> rp_phi(phi_lengths, phi_perm, run_data_phi);
    runperm_phi_inv<RunCols> rp_inv(inv_lengths, inv_perm, run_data_inv);

    using PosPhi = typename runperm_phi<RunCols>::position;
    using PosInv = typename runperm_phi_inv<RunCols>::position;

    PosPhi p0 = rp_phi.first();
    PosPhi phi1 = rp_phi.phi(p0);
    PosPhi next1 = rp_phi.next(p0);
    assert(phi1.interval == next1.interval);
    assert(phi1.offset == next1.offset);

    PosPhi phi3 = rp_phi.phi(p0, 3);
    PosPhi it_phi = p0;
    it_phi = rp_phi.phi(it_phi);
    it_phi = rp_phi.phi(it_phi);
    it_phi = rp_phi.phi(it_phi);
    assert(phi3.interval == it_phi.interval);
    assert(phi3.offset == it_phi.offset);

    PosInv q0 = rp_inv.first();
    PosInv inv1 = rp_inv.phi_inv(q0);
    PosInv next_inv1 = rp_inv.next(q0);
    assert(inv1.interval == next_inv1.interval);
    assert(inv1.offset == next_inv1.offset);

    PosInv inv2 = rp_inv.phi_inv(q0, 2);
    PosInv it_inv = q0;
    it_inv = rp_inv.phi_inv(it_inv);
    it_inv = rp_inv.phi_inv(it_inv);
    assert(inv2.interval == it_inv.interval);
    assert(inv2.offset == it_inv.offset);
}

void test_runperm_phi_inv_img_rank_inv_wrapper_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t inv_domain_1;
    ulint max_length_inv_1;
    auto [inv_lengths_1, inv_perm] = rlbwt_to_phi_inv_images(bwt_heads, bwt_run_lengths, &inv_domain_1, &max_length_inv_1);
    size_t inv_domain_2;
    ulint max_length_inv_2;
    auto [inv_lengths_2, phi_inv_img_rank_inv] = rlbwt_to_phi_inv_img_rank_inv(bwt_heads, bwt_run_lengths, &inv_domain_2, &max_length_inv_2);

    assert(inv_domain_1 == inv_domain_2);
    assert(inv_lengths_1.size() == inv_lengths_2.size());
    assert(inv_perm.size() == phi_inv_img_rank_inv.size());

    for (size_t i = 0; i < inv_lengths_1.size(); ++i) {
        assert(inv_lengths_1[i] == inv_lengths_2[i]);
    }

    auto expected_phi_inv_img_rank_inv = compute_img_rank_inv(inv_perm);
    for (size_t i = 0; i < expected_phi_inv_img_rank_inv.size(); ++i) {
        assert(expected_phi_inv_img_rank_inv[i] == phi_inv_img_rank_inv[i]);
    }
}

void test_runperm_phi_img_rank_inv_wrapper_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain_1;
    ulint max_length_phi_1;
    auto [phi_lengths_1, phi_perm] = rlbwt_to_phi_images(bwt_heads, bwt_run_lengths, &phi_domain_1, &max_length_phi_1);

    size_t phi_domain_2;
    ulint max_length_phi_2;
    auto [phi_lengths_2, phi_img_rank_inv] = rlbwt_to_phi_img_rank_inv(bwt_heads, bwt_run_lengths, &phi_domain_2, &max_length_phi_2);

    assert(phi_domain_1 == phi_domain_2);
    assert(phi_lengths_1.size() == phi_lengths_2.size());
    assert(phi_perm.size() == phi_img_rank_inv.size());

    for (size_t i = 0; i < phi_lengths_1.size(); ++i) {
        assert(phi_lengths_1[i] == phi_lengths_2[i]);
    }

    auto expected_phi_img_rank_inv = compute_img_rank_inv(phi_perm);
    for (size_t i = 0; i < expected_phi_img_rank_inv.size(); ++i) {
        assert(expected_phi_img_rank_inv[i] == phi_img_rank_inv[i]);
    }
}

void test_phi_and_phi_inv_img_rank_inv_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain;
    ulint max_length_phi;
    auto [phi_lengths, phi_img_rank_inv] = rlbwt_to_phi_img_rank_inv(bwt_heads, bwt_run_lengths, &phi_domain, &max_length_phi);
    size_t inv_domain;
    ulint max_length_inv;
    auto [inv_lengths, inv_img_rank_inv] = rlbwt_to_phi_inv_img_rank_inv(bwt_heads, bwt_run_lengths, &inv_domain, &max_length_inv);

    assert(phi_domain == inv_domain);
    assert(phi_lengths.size() == inv_lengths.size());
    assert(phi_img_rank_inv.size() == inv_img_rank_inv.size());

    auto invert_phi_img_rank_inv = get_inverse_permutation(phi_img_rank_inv);
    for (size_t i = 0; i < invert_phi_img_rank_inv.size(); ++i) {
        assert(invert_phi_img_rank_inv[i] == inv_img_rank_inv[i]);
    }
}

int main() {
    test_phi_phi_inv_structure_on_small_rlbwt();
    test_runperm_phi_phi_inv_wrapper_equivalence();
    test_runperm_phi_inv_img_rank_inv_wrapper_equivalence();
    test_runperm_phi_img_rank_inv_wrapper_equivalence();
    test_phi_and_phi_inv_img_rank_inv_equivalence();

    std::cout << "runperm_phi_phi_inv unit tests passed" << std::endl;
    return 0;
}
