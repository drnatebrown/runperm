/* RLBWT benchmarks on real data.
 *
 * Loads precomputed BWT + SA (4-byte little-endian) and runs:
 * - LF / FL text reconstruction (standard + exponential, split + unsplit)
 * - Phi / InvPhi suffix array recovery (standard + exponential, split + unsplit)
 *
 * Datasets (run from the `runperm/` directory):
 * - `tests/data/dna.{txt,bwt,sa}` using `Nucleotide`
 * - `tests/data/hamlet.{txt,bwt,sa}` using generic `Alphabet`
 */

#include "rlbwt.hpp"

#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

using namespace std;
using namespace std::chrono;

// === IO utilities ===

string load_file_as_string(const string &path) {
    ifstream in(path, ios::binary);
    if (!in) {
        throw runtime_error("Failed to open file: " + path);
    }
    string data;
    in.seekg(0, ios::end);
    auto size = in.tellg();
    in.seekg(0, ios::beg);
    data.resize(static_cast<size_t>(size));
    if (size > 0) {
        in.read(&data[0], size);
    }
    return data;
}

vector<uchar> load_file_as_bytes(const string &path) {
    string s = load_file_as_string(path);
    return vector<uchar>(s.begin(), s.end());
}

vector<ulint> load_sa_u32_le(const string &path) {
    ifstream in(path, ios::binary);
    if (!in) {
        throw runtime_error("Failed to open SA file: " + path);
    }
    in.seekg(0, ios::end);
    auto size = in.tellg();
    in.seekg(0, ios::beg);
    if (size < 0 || (static_cast<size_t>(size) % 4) != 0) {
        throw runtime_error("SA file size is not a multiple of 4 bytes: " + path);
    }
    size_t n = static_cast<size_t>(size) / 4;
    vector<ulint> sa;
    sa.resize(n);
    for (size_t i = 0; i < n; ++i) {
        uint32_t v = 0;
        in.read(reinterpret_cast<char*>(&v), sizeof(v));
        if (!in) {
            throw runtime_error("Failed while reading SA file: " + path);
        }
        sa[i] = static_cast<ulint>(v);
    }
    return sa;
}

void ensure_nucleotide_bytes(const vector<uchar>& bytes, const string& what) {
    auto ok = [](uchar c) {
        return c == 0 || c == 1 || c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
    };
    vector<uchar> bad;
    bad.reserve(8);
    for (uchar c : bytes) {
        if (!ok(c)) bad.push_back(c);
    }
    if (!bad.empty()) {
        sort(bad.begin(), bad.end());
        bad.erase(unique(bad.begin(), bad.end()), bad.end());
        ostringstream oss;
        oss << "Non-nucleotide byte(s) found in " << what << ":";
        for (uchar c : bad) {
            oss << " 0x" << hex << static_cast<int>(c);
        }
        throw runtime_error(oss.str());
    }
}

// Build a simple suffix array using the classic doubling algorithm.
// This is only for benchmarking correctness vs the RLBWT-based SA.
vector<ulint> build_suffix_array(const string &s) {
    const size_t n = s.size();
    vector<int> sa(n);
    vector<int> rank(n + 1);
    vector<int> tmp(n + 1);

    for (size_t i = 0; i < n; ++i) {
        sa[i] = static_cast<int>(i);
        rank[i] = static_cast<unsigned char>(s[i]);
    }
    rank[n] = -1;

    for (size_t k = 1; k < n; k <<= 1) {
        auto cmp = [&](int i, int j) {
            if (rank[i] != rank[j]) return rank[i] < rank[j];
            int ri = (i + static_cast<int>(k) <= static_cast<int>(n)) ? rank[i + static_cast<int>(k)] : -1;
            int rj = (j + static_cast<int>(k) <= static_cast<int>(n)) ? rank[j + static_cast<int>(k)] : -1;
            return ri < rj;
        };
        sort(sa.begin(), sa.end(), cmp);
        tmp[sa[0]] = 0;
        for (size_t i = 1; i < n; ++i) {
            tmp[sa[i]] = tmp[sa[i - 1]] + (cmp(sa[i - 1], sa[i]) ? 1 : 0);
        }
        for (size_t i = 0; i < n; ++i) {
            rank[i] = tmp[i];
        }
        if (rank[sa[n - 1]] == static_cast<int>(n - 1)) break;
    }

    vector<ulint> sa_ul(n);
    for (size_t i = 0; i < n; ++i) {
        sa_ul[i] = static_cast<ulint>(sa[i]);
    }
    return sa_ul;
}

// === LF / FL text reconstruction benchmarks ===

template<typename MoveLFType>
void bench_move_lf(const string &name,
                   const vector<uchar> &bwt_heads,
                   const vector<ulint> &bwt_run_lengths,
                   const string &text,
                   const SplitParams &split_params) {
    cout << "  " << name << " (LF)" << endl;

    auto t0 = high_resolution_clock::now();
    MoveLFType move_lf(bwt_heads, bwt_run_lengths, split_params);
    auto t1 = high_resolution_clock::now();

    using Position = typename MoveLFType::Position;
    string recovered_text(text.size(), '\0');
    Position pos = move_lf.first();
    for (size_t i = 1; i < move_lf.domain(); ++i) {
        recovered_text[text.size() - i] = static_cast<char>(move_lf.get_character(pos));
        pos = move_lf.LF(pos);
    }
    assert(recovered_text.compare(text) == 0);
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto reconstruction_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    Reconstruction: " << reconstruction_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_lf.serialize(ss) << " bytes" << endl;
}

template<typename MoveFLType>
void bench_move_fl(const string &name,
                   const vector<uchar> &bwt_heads,
                   const vector<ulint> &bwt_run_lengths,
                   const string &text,
                   const SplitParams &split_params) {
    cout << "  " << name << " (FL)" << endl;

    auto t0 = high_resolution_clock::now();
    MoveFLType move_fl(bwt_heads, bwt_run_lengths, split_params);
    auto t1 = high_resolution_clock::now();

    using Position = typename MoveFLType::Position;
    string recovered_text = "";
    Position pos = move_fl.first();
    pos = move_fl.FL(pos);
    for (size_t i = 1; i < move_fl.domain(); ++i) {
        recovered_text += static_cast<char>(move_fl.get_character(pos));
        pos = move_fl.FL(pos);
    }
    assert(recovered_text.compare(text) == 0);
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto reconstruction_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    Reconstruction: " << reconstruction_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_fl.serialize(ss) << " bytes" << endl;
}

template<typename MoveLFStd, typename MoveLFExp, typename MoveFLStd, typename MoveFLExp>
void run_lf_fl_benchmarks(const vector<uchar> &bwt_heads,
                          const vector<ulint> &bwt_run_lengths,
                          const string &text) {
    cout << "=== LF / FL reconstruction benchmarks ===" << endl;
    cout << "Text length: " << text.size() << endl;
    cout << "RLBWT runs:  " << bwt_heads.size() << endl;

    SplitParams no_splitting = NO_SPLITTING;
    SplitParams splitting(8.0, std::nullopt); // length_capping_factor = 8.0

    // LF (standard + exponential)
    bench_move_lf<MoveLFStd>("MoveLF (no splitting)", bwt_heads, bwt_run_lengths, text, no_splitting);
    bench_move_lf<MoveLFStd>("MoveLF (split, 8.0)", bwt_heads, bwt_run_lengths, text, splitting);
    bench_move_lf<MoveLFExp>("MoveLFExp (no splitting)", bwt_heads, bwt_run_lengths, text, no_splitting);
    bench_move_lf<MoveLFExp>("MoveLFExp (split, 8.0)", bwt_heads, bwt_run_lengths, text, splitting);

    // FL (standard + exponential)
    bench_move_fl<MoveFLStd>("MoveFL (no splitting)", bwt_heads, bwt_run_lengths, text, no_splitting);
    bench_move_fl<MoveFLStd>("MoveFL (split, 8.0)", bwt_heads, bwt_run_lengths, text, splitting);
    bench_move_fl<MoveFLExp>("MoveFLExp (no splitting)", bwt_heads, bwt_run_lengths, text, no_splitting);
    bench_move_fl<MoveFLExp>("MoveFLExp (split, 8.0)", bwt_heads, bwt_run_lengths, text, splitting);

    cout << endl;
}

// === Phi / InvPhi suffix array benchmarks (move-only, mirroring tests/integration/rlbwt_test.cpp) ===

// Standard Phi (no exponential search)
template<class Container1, class Container2>
void bench_move_phi(const string &name,
                    const Container1 &lengths,
                    const Container2 &tau_inv,
                    ulint domain,
                    const vector<ulint> &sa_truth,
                    const SplitParams &split_params) {
    cout << "  " << name << endl;

    auto t0 = high_resolution_clock::now();
    auto permutation = Permutation::from_lengths_and_tau_inv(lengths, tau_inv);
    MovePhi move_phi(permutation);
    auto t1 = high_resolution_clock::now();

    using Position = typename MovePhi::Position;
    vector<ulint> sa_recovered(sa_truth.size());

    Position pos = move_phi.last();
    pos = move_phi.Phi(pos);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[sa_recovered.size() - i - 1] = move_phi.SA(pos);
        pos = move_phi.Phi(pos);
    }

    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa_truth[i]);
    }
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto sa_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    SA recovery:    " << sa_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_phi.serialize(ss) << " bytes" << endl;
}

// Phi with exponential search (absolute positions + ExponentialSearch = true)
template<class Container1, class Container2>
void bench_move_phi_exp(const string &name,
                        const Container1 &lengths,
                        const Container2 &tau_inv,
                        ulint domain,
                        const vector<ulint> &sa_truth,
                        const SplitParams &split_params) {
    cout << "  " << name << endl;

    auto t0 = high_resolution_clock::now();
    class MovePhiExp : public MovePermImpl<true, true, MoveCols, MoveStructure, MoveVector> {
        using Base = MovePermImpl<true, true, MoveCols, MoveStructure, MoveVector>;
    public:
        using Base::Base;
        using Base::operator=;
        using Position = typename Base::Position;

        Position Phi(Position pos) { return Base::next(pos); }
        Position Phi(Position pos, ulint steps) { return Base::next(pos, steps); }
        ulint SA(Position pos) { return pos.idx; }
    };

    auto permutation = Permutation::from_lengths_and_tau_inv(lengths, tau_inv);
    MovePhiExp move_phi(permutation);
    auto t1 = high_resolution_clock::now();

    using Position = typename MovePhiExp::Position;
    vector<ulint> sa_recovered(sa_truth.size());

    Position pos = move_phi.last();
    pos = move_phi.Phi(pos);
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[sa_recovered.size() - i - 1] = move_phi.SA(pos);
        pos = move_phi.Phi(pos);
    }

    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa_truth[i]);
    }
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto sa_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    SA recovery:    " << sa_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_phi.serialize(ss) << " bytes" << endl;
}

// Standard InvPhi (no exponential search)
template<class Container1, class Container2>
void bench_move_invphi(const string &name,
                       const Container1 &lengths,
                       const Container2 &tau_inv,
                       ulint domain,
                       const vector<ulint> &sa_truth,
                       const SplitParams &split_params) {
    cout << "  " << name << endl;

    auto t0 = high_resolution_clock::now();
    auto permutation = Permutation::from_lengths_and_tau_inv(lengths, tau_inv);
    MoveInvPhi move_invphi(permutation);
    auto t1 = high_resolution_clock::now();

    using Position = typename MoveInvPhi::Position;
    vector<ulint> sa_recovered(sa_truth.size());

    Position pos = move_invphi.last();
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[i] = move_invphi.SA(pos);
        pos = move_invphi.InvPhi(pos);
    }

    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa_truth[i]);
    }
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto sa_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    SA recovery:    " << sa_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_invphi.serialize(ss) << " bytes" << endl;
}

// InvPhi with exponential search (absolute positions + ExponentialSearch = true)
template<class Container1, class Container2>
void bench_move_invphi_exp(const string &name,
                           const Container1 &lengths,
                           const Container2 &tau_inv,
                           ulint domain,
                           const vector<ulint> &sa_truth,
                           const SplitParams &split_params) {
    cout << "  " << name << endl;

    auto t0 = high_resolution_clock::now();
    class MoveInvPhiExp : public MovePermImpl<true, true, MoveCols, MoveStructure, MoveVector> {
        using Base = MovePermImpl<true, true, MoveCols, MoveStructure, MoveVector>;
    public:
        using Base::Base;
        using Base::operator=;
        using Position = typename Base::Position;

        Position InvPhi(Position pos) { return Base::next(pos); }
        Position InvPhi(Position pos, ulint steps) { return Base::next(pos, steps); }
        ulint SA(Position pos) { return pos.idx; }
    };

    auto permutation = Permutation::from_lengths_and_tau_inv(lengths, tau_inv);
    MoveInvPhiExp move_invphi(permutation);
    auto t1 = high_resolution_clock::now();

    using Position = typename MoveInvPhiExp::Position;
    vector<ulint> sa_recovered(sa_truth.size());

    Position pos = move_invphi.last();
    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        sa_recovered[i] = move_invphi.SA(pos);
        pos = move_invphi.InvPhi(pos);
    }

    for (size_t i = 0; i < sa_recovered.size(); ++i) {
        assert(sa_recovered[i] == sa_truth[i]);
    }
    auto t3 = high_resolution_clock::now();

    auto creation_duration = duration_cast<microseconds>(t1 - t0);
    auto sa_duration = duration_cast<microseconds>(t3 - t1);
    auto total_duration = duration_cast<microseconds>(t3 - t0);

    cout << "    Creation:       " << creation_duration.count() << "us" << endl;
    cout << "    SA recovery:    " << sa_duration.count() << "us" << endl;
    cout << "    Total:          " << total_duration.count() << "us" << endl;
    stringstream ss;
    cout << "    Size:           " << move_invphi.serialize(ss) << " bytes" << endl;
}

template<typename AlphabetType=Nucleotide>
void run_phi_invphi_benchmarks(const vector<uchar> &bwt_heads,
                               const vector<ulint> &bwt_run_lengths,
                               const vector<ulint> &sa_truth) {
    cout << "=== Phi / InvPhi SA benchmarks ===" << endl;

    cout << "Domain (SA size): " << sa_truth.size() << endl;

    // Build Phi / InvPhi structures from the RLBWT.
    size_t phi_domain;
    ulint max_length;
    auto t_phi_derive = high_resolution_clock::now();
    auto [phi_lengths, phi_tau_inv] = phi::rlbwt_to_phi_tau_inv<AlphabetType>(bwt_heads, bwt_run_lengths, &phi_domain, &max_length);
    auto t_phi_derive_end = high_resolution_clock::now();
    auto phi_derive_duration = duration_cast<microseconds>(t_phi_derive_end - t_phi_derive);
    cout << "    Phi derivation: " << phi_derive_duration.count() << "us" << endl;

    size_t invphi_domain;
    ulint max_length_inv;
    auto t_invphi_derive = high_resolution_clock::now();
    auto [invphi_lengths, invphi_tau_inv] = invphi::rlbwt_to_invphi_tau_inv<AlphabetType>(bwt_heads, bwt_run_lengths, &invphi_domain, &max_length_inv);
    auto t_invphi_derive_end = high_resolution_clock::now();
    auto invphi_derive_duration = duration_cast<microseconds>(t_invphi_derive_end - t_invphi_derive);
    cout << "    InvPhi derivation: " << invphi_derive_duration.count() << "us" << endl;

    assert(phi_domain == sa_truth.size());
    assert(invphi_domain == sa_truth.size());

    SplitParams no_splitting = NO_SPLITTING;
    SplitParams splitting(8.0, std::nullopt); // length_capping_factor = 8.0

    // Move-only Phi (standard + exponential)
    bench_move_phi(
        "MovePhi (no splitting)",
        phi_lengths,
        phi_tau_inv,
        phi_domain,
        sa_truth,
        no_splitting);

    bench_move_phi(
        "MovePhi (split, 8.0)",
        phi_lengths,
        phi_tau_inv,
        phi_domain,
        sa_truth,
        splitting);
    bench_move_phi_exp(
        "MovePhiExp (no splitting)",
        phi_lengths,
        phi_tau_inv,
        phi_domain,
        sa_truth,
        no_splitting);
    bench_move_phi_exp(
        "MovePhiExp (split, 8.0)",
        phi_lengths,
        phi_tau_inv,
        phi_domain,
        sa_truth,
        splitting);

    // Move-only InvPhi (standard + exponential)
    bench_move_invphi(
        "MoveInvPhi (no splitting)",
        invphi_lengths,
        invphi_tau_inv,
        invphi_domain,
        sa_truth,
        no_splitting);

    bench_move_invphi(
        "MoveInvPhi (split, 8.0)",
        invphi_lengths,
        invphi_tau_inv,
        invphi_domain,
        sa_truth,
        splitting);

    bench_move_invphi_exp(
        "MoveInvPhiExp (no splitting)",
        invphi_lengths,
        invphi_tau_inv,
        invphi_domain,
        sa_truth,
        no_splitting);

    bench_move_invphi_exp(
        "MoveInvPhiExp (split, 8.0)",
        invphi_lengths,
        invphi_tau_inv,
        invphi_domain,
        sa_truth,
        splitting);

    cout << endl;
}

int main() {
    try {
        // Paths are resolved relative to the process working directory.
        // Run this benchmark from the `runperm/` directory.
        struct Dataset {
            string name;
            string text_path;
            string bwt_path;
            string sa_path;
            bool nucleotide;
        };

        vector<Dataset> datasets = {
            {"dna (nucleotide)",   "tests/data/dna.txt",    "tests/data/dna.bwt",    "tests/data/dna.sa",    true},
            {"hamlet (alphabet)",  "tests/data/hamlet.txt", "tests/data/hamlet.bwt", "tests/data/hamlet.sa", false},
        };

        for (const auto& ds : datasets) {
            cout << "==============================" << endl;
            cout << "Dataset: " << ds.name << endl;
            cout << "Text:    " << ds.text_path << endl;
            cout << "BWT:     " << ds.bwt_path << endl;
            cout << "SA:      " << ds.sa_path << endl;

            string text = load_file_as_string(ds.text_path);
            vector<uchar> bwt_chars = load_file_as_bytes(ds.bwt_path);
            vector<ulint> sa_truth = load_sa_u32_le(ds.sa_path);

            if (bwt_chars.size() != sa_truth.size()) {
                throw runtime_error("BWT length != SA length for dataset: " + ds.name);
            }
            if (bwt_chars.size() != text.size() + 1) {
                throw runtime_error("Expected |BWT| == |text| + 1 (NUL terminator) for dataset: " + ds.name);
            }

            if (ds.nucleotide) {
                ensure_nucleotide_bytes(bwt_chars, ds.bwt_path);
                ensure_nucleotide_bytes(vector<uchar>(text.begin(), text.end()), ds.text_path);
            }

            auto [bwt_heads, bwt_run_lengths] = bwt_to_rlbwt(bwt_chars);

            cout << "Loaded text length: " << text.size() << endl;
            cout << "Loaded BWT length:  " << bwt_chars.size() << endl;
            cout << "RLBWT runs:         " << bwt_heads.size() << endl;
            cout << endl;

            if (ds.nucleotide) {
                using LFStd = MoveLF<>;
                using LFExp = MoveLFImpl<true, true>;
                using FLStd = MoveFL<>;
                using FLExp = MoveFLImpl<true, true>;
                run_lf_fl_benchmarks<LFStd, LFExp, FLStd, FLExp>(bwt_heads, bwt_run_lengths, text);

                using LFPhi = MoveLF<>;
                run_phi_invphi_benchmarks<Nucleotide>(bwt_heads, bwt_run_lengths, sa_truth);
            } else {
                using LFStd = MoveLFImpl<false, false, Alphabet>;
                using LFExp = MoveLFImpl<true, true, Alphabet>;
                using FLStd = MoveFLImpl<false, false, Alphabet>;
                using FLExp = MoveFLImpl<true, true, Alphabet>;
                run_lf_fl_benchmarks<LFStd, LFExp, FLStd, FLExp>(bwt_heads, bwt_run_lengths, text);

                using LFPhi = MoveLFImpl<false, false, Alphabet>;
                run_phi_invphi_benchmarks<Alphabet>(bwt_heads, bwt_run_lengths, sa_truth);
            }
        }

        cout << "All datasets complete." << endl;
        return 0;
    } catch (const exception &ex) {
        cerr << "Error in rlbwt_bench: " << ex.what() << endl;
        return 1;
    }
}
