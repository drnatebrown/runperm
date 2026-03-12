// Tiny unit tests for Phi / InvPhi helpers and wrappers.
// Focus: structural sanity and wrapper behaviour, not full SA reconstruction.

#include "rlbwt.hpp"

#include <cassert>
#include <iostream>
#include <numeric>
#include <vector>

using std::size_t;
using std::vector;

void test_phi_invphi_structure_on_small_rlbwt() {
    // Use the same small example as in the integration test, which is known to be valid.
    vector<uchar> heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> lens  =       { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain;
    ulint max_length;
    auto [phi_lengths, phi_perm] = phi::rlbwt_to_phi(heads, lens, &phi_domain, &max_length);
    size_t inv_domain;
    ulint max_length_inv;
    auto [inv_lengths, inv_perm] = phi::rlbwt_to_invphi(heads, lens, &inv_domain, &max_length_inv);

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

void test_runperm_phi_invphi_wrapper_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t phi_domain;
    ulint max_length;
    auto [phi_lengths, phi_perm] = phi::rlbwt_to_phi(bwt_heads, bwt_run_lengths, &phi_domain, &max_length);
    size_t inv_domain;
    ulint max_length_inv;
    auto [inv_lengths, inv_perm] = phi::rlbwt_to_invphi(bwt_heads, bwt_run_lengths, &inv_domain, &max_length_inv);

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

    RunPermPhi<RunCols> rp_phi(phi_lengths, phi_perm, run_data_phi);
    RunPermInvPhi<RunCols> rp_inv(inv_lengths, inv_perm, run_data_inv);

    using PosPhi = typename RunPermPhi<RunCols>::Position;
    using PosInv = typename RunPermInvPhi<RunCols>::Position;

    PosPhi p0 = rp_phi.first();
    PosPhi phi1 = rp_phi.Phi(p0);
    PosPhi next1 = rp_phi.next(p0);
    assert(phi1.interval == next1.interval);
    assert(phi1.offset == next1.offset);

    PosPhi phi3 = rp_phi.Phi(p0, 3);
    PosPhi it_phi = p0;
    it_phi = rp_phi.Phi(it_phi);
    it_phi = rp_phi.Phi(it_phi);
    it_phi = rp_phi.Phi(it_phi);
    assert(phi3.interval == it_phi.interval);
    assert(phi3.offset == it_phi.offset);

    PosInv q0 = rp_inv.first();
    PosInv inv1 = rp_inv.InvPhi(q0);
    PosInv next_inv1 = rp_inv.next(q0);
    assert(inv1.interval == next_inv1.interval);
    assert(inv1.offset == next_inv1.offset);

    PosInv inv2 = rp_inv.InvPhi(q0, 2);
    PosInv it_inv = q0;
    it_inv = rp_inv.InvPhi(it_inv);
    it_inv = rp_inv.InvPhi(it_inv);
    assert(inv2.interval == it_inv.interval);
    assert(inv2.offset == it_inv.offset);
}

void test_runperm_invphi_tau_inv_wrapper_equivalence() {
    // Reuse the small RLBWT example from the integration test.
    vector<uchar> bwt_heads =       {'T','C','G','A','T', 1 ,'A','T','A'};
    vector<ulint> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

    size_t inv_domain_1;
    ulint max_length_inv_1;
    auto [inv_lengths_1, inv_perm] = phi::rlbwt_to_invphi(bwt_heads, bwt_run_lengths, &inv_domain_1, &max_length_inv_1);
    size_t inv_domain_2;
    ulint max_length_inv_2;
    auto [inv_lengths_2, invphi_tau_inv] = phi::rlbwt_to_invphi_tau_inv(bwt_heads, bwt_run_lengths, &inv_domain_2, &max_length_inv_2);

    assert(inv_domain_1 == inv_domain_2);
    assert(inv_lengths_1.size() == inv_lengths_2.size());
    assert(inv_perm.size() == invphi_tau_inv.size());

    for (size_t i = 0; i < inv_lengths_1.size(); ++i) {
        assert(inv_lengths_1[i] == inv_lengths_2[i]);
    }

    auto expected_invphi_tau_inv = compute_tau_inv(inv_perm);
    for (size_t i = 0; i < expected_invphi_tau_inv.size(); ++i) {
        assert(expected_invphi_tau_inv[i] == invphi_tau_inv[i]);
    }
}

int main() {
    test_phi_invphi_structure_on_small_rlbwt();
    test_runperm_phi_invphi_wrapper_equivalence();
    test_runperm_invphi_tau_inv_wrapper_equivalence();

    std::cout << "runperm_phi_invphi unit tests passed" << std::endl;
    return 0;
}
