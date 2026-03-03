// Unit tests for Alphabet and Nucleotide helper.
// These are simple assert-based tests, no external framework.

#include "internal/ds/alphabet.hpp"

#include <cassert>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::size_t;
using std::vector;

void test_alphabet_basic_mapping() {
    vector<ulint> char_counts(MAX_ALPHABET_SIZE, 0);

    // Enable a small subset of characters.
    char_counts['A'] = 3;
    char_counts['C'] = 5;
    char_counts['Z'] = 1;

    Alphabet alphabet(char_counts);

    assert(alphabet.size() == 3);

    uchar a_code = alphabet.map_char(static_cast<uchar>('A'));
    uchar c_code = alphabet.map_char(static_cast<uchar>('C'));
    uchar z_code = alphabet.map_char(static_cast<uchar>('Z'));

    // Mapping is assigned in increasing character order.
    assert(a_code == 0);
    assert(c_code == 1);
    assert(z_code == 2);

    assert(alphabet.unmap_char(a_code) == static_cast<uchar>('A'));
    assert(alphabet.unmap_char(c_code) == static_cast<uchar>('C'));
    assert(alphabet.unmap_char(z_code) == static_cast<uchar>('Z'));

    vector<uchar> seq{
        static_cast<uchar>('A'),
        static_cast<uchar>('Z'),
        static_cast<uchar>('C')
    };
    auto mapped = alphabet.map_sequence(seq);
    auto unmapped = alphabet.unmap_sequence(mapped);

    assert(unmapped.size() == seq.size());
    for (size_t i = 0; i < seq.size(); ++i) {
        assert(unmapped[i] == seq[i]);
    }
}

void test_alphabet_serialize_roundtrip() {
    vector<ulint> char_counts(MAX_ALPHABET_SIZE, 0);
    char_counts['A'] = 1;
    char_counts['B'] = 2;
    char_counts['C'] = 3;

    Alphabet alphabet(char_counts);

    std::stringstream ss;
    size_t bytes_written = alphabet.serialize(ss);
    assert(bytes_written > 0);

    Alphabet loaded;
    loaded.load(ss);

    // Check that size and mappings for the active characters are preserved.
    assert(loaded.size() == alphabet.size());

    uchar a_code = loaded.map_char(static_cast<uchar>('A'));
    uchar b_code = loaded.map_char(static_cast<uchar>('B'));
    uchar c_code = loaded.map_char(static_cast<uchar>('C'));

    assert(loaded.unmap_char(a_code) == static_cast<uchar>('A'));
    assert(loaded.unmap_char(b_code) == static_cast<uchar>('B'));
    assert(loaded.unmap_char(c_code) == static_cast<uchar>('C'));
}

void test_alphabet_full_alphabet_throws() {
    vector<ulint> char_counts(MAX_ALPHABET_SIZE, 1);

    bool threw = false;
    try {
        Alphabet alphabet(char_counts);
    } catch (const std::runtime_error &) {
        threw = true;
    }
    assert(threw);
}

void test_nucleotide_basic_mapping() {
    vector<ulint> char_counts(MAX_ALPHABET_SIZE, 0);

    // Terminator and separator plus the standard DNA alphabet.
    char_counts[TERMINATOR] = 1;
    char_counts[SEPARATOR] = 1;
    char_counts[static_cast<uchar>('A')] = 1;
    char_counts[static_cast<uchar>('C')] = 1;
    char_counts[static_cast<uchar>('G')] = 1;
    char_counts[static_cast<uchar>('T')] = 1;
    char_counts[static_cast<uchar>('N')] = 1;

    Nucleotide nuc(char_counts);
    (void)nuc; // suppress unused warning in release builds

    assert(Nucleotide::size() == Nucleotide::SIGMA);

    // Terminator and separator map to the first two codes.
    assert(Nucleotide::map_char(static_cast<uchar>(TERMINATOR)) == 0);
    assert(Nucleotide::map_char(static_cast<uchar>(SEPARATOR)) == 1);

    // Upper- and lower-case DNA symbols should map consistently.
    uchar a_code = Nucleotide::map_char(static_cast<uchar>('A'));
    uchar c_code = Nucleotide::map_char(static_cast<uchar>('C'));
    uchar g_code = Nucleotide::map_char(static_cast<uchar>('G'));
    uchar t_code = Nucleotide::map_char(static_cast<uchar>('T'));
    uchar n_code = Nucleotide::map_char(static_cast<uchar>('N'));

    uchar a_lower_code = Nucleotide::map_char(static_cast<uchar>('a'));
    uchar c_lower_code = Nucleotide::map_char(static_cast<uchar>('c'));
    uchar g_lower_code = Nucleotide::map_char(static_cast<uchar>('g'));
    uchar t_lower_code = Nucleotide::map_char(static_cast<uchar>('t'));
    uchar n_lower_code = Nucleotide::map_char(static_cast<uchar>('n'));

    assert(a_code == a_lower_code);
    assert(c_code == c_lower_code);
    assert(g_code == g_lower_code);
    assert(t_code == t_lower_code);
    assert(n_code == n_lower_code);

    // Reverse mapping yields canonical (upper-case) characters.
    assert(Nucleotide::unmap_char(a_code) == static_cast<uchar>('A'));
    assert(Nucleotide::unmap_char(c_code) == static_cast<uchar>('C'));
    assert(Nucleotide::unmap_char(g_code) == static_cast<uchar>('G'));
    assert(Nucleotide::unmap_char(t_code) == static_cast<uchar>('T'));
    assert(Nucleotide::unmap_char(n_code) == static_cast<uchar>('N'));
}

void test_nucleotide_sequence_roundtrip() {
    vector<ulint> char_counts(MAX_ALPHABET_SIZE, 0);
    char_counts[TERMINATOR] = 1;
    char_counts[SEPARATOR] = 1;
    char_counts[static_cast<uchar>('A')] = 1;
    char_counts[static_cast<uchar>('C')] = 1;
    char_counts[static_cast<uchar>('G')] = 1;
    char_counts[static_cast<uchar>('T')] = 1;
    char_counts[static_cast<uchar>('N')] = 1;

    Nucleotide nuc(char_counts);
    (void)nuc;

    vector<uchar> seq{
        static_cast<uchar>('A'),
        static_cast<uchar>('c'),
        static_cast<uchar>('G'),
        static_cast<uchar>('t'),
        static_cast<uchar>('N')
    };

    auto mapped = Nucleotide::map_sequence(seq);
    auto unmapped = Nucleotide::unmap_sequence(mapped);

    // Lowercase inputs should be normalized to uppercase in the unmapped sequence.
    vector<uchar> expected{
        static_cast<uchar>('A'),
        static_cast<uchar>('C'),
        static_cast<uchar>('G'),
        static_cast<uchar>('T'),
        static_cast<uchar>('N')
    };

    assert(unmapped.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        assert(unmapped[i] == expected[i]);
    }
}

int main() {
    test_alphabet_basic_mapping();
    test_alphabet_serialize_roundtrip();
    test_alphabet_full_alphabet_throws();

    test_nucleotide_basic_mapping();
    test_nucleotide_sequence_roundtrip();

    std::cout << "alphabet tests passed" << std::endl;
    return 0;
}

