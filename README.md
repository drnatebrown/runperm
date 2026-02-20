# RunPerm

Flexible plug-and-play header library implementing data structures for representing run-length encoded permutations. Features move structures and specialized RLBWT implementations for high-performance indexing. Can be used as a foundation for further applications, providing storage and retrieval of user fields alongside the permutation intervals.

Described in greater detail at [https://arxiv.org/abs/2602.11029](https://arxiv.org/abs/2602.11029).

Alpha Version 0.1.3

## Overview

RunPerm provides compact representations of run-length encoded permutations, with specialized implementations for various indexing applications. The core permutation logic is stored in a move structure. The library offers three main interfaces:

- **RunPerm**: Core run-length permutation data structure allowing user defined fields stored alongside permutation intervals.
- **MovePerm**: Simplified permutation representation using move structures, no user data.
- **RLBWT**: Specialized run-length Burrows-Wheeler Transform permutations (LF/FL and Phi/InvPhi) built on RunPerm and MovePerm.

## Features

- **Compact Representation**: Bitpacked for minimum fixed width per stored value.
- **Length capping**: Optimization scheme for faster and smaller move structures.
- **Move Structures**: Advanced data structures for fast and cache efficient navigation of permutations.
- **RLBWT Specializations**:
  - LF/FL navigation with character access.
  - Phi/InvPhi navigation with SA retrieval.
  - Helpers to derive Phi/InvPhi structures from run-length BWT.
- **Flexible Configuration**: Multiple template parameters for various move structure representations.
- **Advanced Storage** Allows efficient storage and retrieval of additional information alongside permutations.

## Quick Start

### Installation

This library is header-only. Just add `include/` to your compiler's include paths and include `runperm.hpp`, `move.hpp`, or `rlbwt.hpp` as needed:
- `runperm.hpp` loads RunPerm and MovePerm.
- `move.hpp` loads the underlying move structure implementation only.
- `rlbwt.hpp` loads specialized classes to represent LF/FL and Phi/InvPhi.

The provided `makefile` only builds example/test executables from `src`:
- `example_test`
- `move_test`
- `rlbwt_test`
- `runperm_test`

### Basic Usage

#### RunPerm
`RunPerm` is the generic implementation allowing storage of user defined fields alongside a runny permutation, included with `runperm.hpp`. Some small examples are given in the `src/example_test.cpp` file.

Given the runny permutation example:

<img width="598" height="572" alt="Image" src="https://github.com/user-attachments/assets/18933a68-b858-40eb-bfb3-f5fd758557d6" />

We pass the lengths of contigiously permuted intervals and the permutation of their first values:

```cpp
#include "runperm.hpp"

// Create a run-length permutation
// Length of contigiuously permuted intervals
std::vector<unsigned long int> lengths = {2, 3, 1, 2, 2, 1, 1, 1, 3};
// Permutation at head of each intervals
std::vector<unsigned long int> permutation = {1, 9, 3, 12, 4, 14, 0, 15, 6};   
ulint domain = 16;                                // Total domain size

// Some example data to store alongside these runs.
// Must always include COUNT as the last entry to signal number of fields.
enum class RunCols {
    VAL1,
    VAL2,
    COUNT
};
// Defines an array of unsigned long integers of length COUNT
using RunColsTuple = DataTuple<RunCols>;
std::vector<RunColsTuple> run_data(lengths.size()); // Some data with tuples per run

// Basic construction
RunPerm<RunCols> rp(lengths, permutation, domain, run_data);
```

#### MovePerm
`MovePerm` functions similarly to `RunPerm` but without user defined fields, also included in `runperm.hpp`

```cpp
#include "runperm.hpp"

// Create from permutation vector
std::vector<unsigned long int> permutation = {1, 2, 9, 10, 11, 3, 12, 13, 4, 5, 14, 0, 15, 6, 7, 8};
// Stores interval/offset pairs as positions
MovePermRelative mp(permutation); 

// Or from run structure
auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
ulint domain = permutation.size();
// Also stores the absolute position in the permutation
MovePermAbsolute mp_abs(lengths, permutation, domain);
```

#### Length Capping
TODO

#### RLBWT: LF/FL
By including `rlbwt.hpp` users can use specialized methods designed for permutations based on the RLBWT such as LF/FL for pattern matching and text extraction.

```cpp
#include "rlbwt.hpp"

// BWT = TTTTTCCCGGGAAAT$ATTTTAAAAAA
// RLBWT as run heads (characters) and run lengths, with 0 as terminator
std::vector<uchar> bwt_heads =       {'T','C','G','A','T', 0 ,'A','T','A'};
std::vector<unsigned long int> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

// LF using Move-only interface
MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

// LF with RunPerm + run data columns
enum class RunCols { VAL1, VAL2, COUNT };
using LFRunData = DataTuple<RunCols>;
std::vector<LFRunData> lf_run_data(bwt_heads.size()); // Insert with some type of data
RunPermLF<RunCols> runperm_lf(bwt_heads, bwt_run_lengths, lf_run_data);

// FL (symmetrically)
MoveFL<> move_fl(bwt_heads, bwt_run_lengths);
...
```

#### RLBWT: Phi/InvPhi
By including `rlbwt.hpp`, we can also build the permutations from an RLBWT for $\phi$ and $\phi^{-1}$ needed for locate queries.

```cpp
#include "rlbwt.hpp"

// Build Phi structure from RLBWT
auto [phi_lengths, phi_interval_perm, domain] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);

// RunPerm Phi (always absolute positions)
enum class PhiCols { VAL1, VAL2, COUNT };
using PhiRunData = DataTuple<RunCols>;
std::vector<PhiRunData> phi_run_data(phi_lengths.size());
RunPermPhi<PhiCols> phi(phi_lengths, phi_interval_perm, domain, phi_run_data);

// Move-only Phi
auto [phi_lens2, phi_perm2, domain2] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
MovePhi mphi(phi_lens2, phi_perm2, domain2);

// Build InvPhi similarly...
auto [invphi_lengths, invphi_interval_perm, domain_inv] = rlbwt_to_invphi(bwt_heads, bwt_run_lengths);
...

```

## Interface

### RunPerm

The main class for run-length encoded permutations with move structure integration.

- **Template Parameters**:
  - `RunColsType`: Type defining user run data columns
  - `IntegratedMoveStructure`: integrate run data with move structure (default: false)
  - `StoreAbsolutePositions`: store absolute positions for index lookups (default: false)
- **Key Methods**:
  - `next(pos)`, `next(pos, ulint steps)`
  - `up(pos)`, `down(pos)`, `first()`, `last()`
  - `get<Col>(pos)`, `get<Col>(pos, i)`, `get_length(pos)`, `get_length(i)`
  - `pred<Col>(pos, val)`, `succ<Col>(pos, val)`
  - `domain()`, `move_runs()`, `permutation_runs()`
  - `serialize(os)`, `load(is)`

### Performance Considerations
- **Integrated vs Separated**: Integrating user data alongside the move structure offer better cache locality but may cause slower move queries since navigating the data structure requires loading larger entries.
- **Absolute vs Relative Positions**: Absolute positions enable full permutation positional information but increase memory usage. The space usage is, for a runny permutation of $r$ runs over domain $n$ is approximately:
  - **Absolute**: $r \log r + 2 r \log n$ bits
  - **Relative**: $r \log r + 2 r \log \frac{n}{r}$ bits

### Dependencies

- C++17 or later
- No external dependencies; header-only usage.

## Examples

### Basic Permutation Navigation

```cpp
#include "runperm.hpp"

RunPerm<MyRunCols> rp(lengths, permutation, domain, run_data);

// Navigate
auto pos = rp.first();
pos = rp.next(pos);
pos = rp.next(pos, 5);
pos = rp.down(pos);
```

### RLBWT LF Text Recovery

```cpp
#include "rlbwt.hpp"

MoveLF<> lf(bwt_heads, bwt_run_lengths);
pos = lf.first();
std::string recovered(text.size(), '\0');
for (size_t i = 1; i < lf.domain(); ++i) {
    recovered[text.size() - i] = (char)lf.get_character(pos);
    pos = lf.LF(pos);
}
```

### Phi SA Recovery

```cpp
#include "rlbwt.hpp"

auto [phi_lens, phi_perm, n] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
MovePhi phi(phi_lens, phi_perm, n);
auto pos = phi.last();
pos = phi.Phi(pos);
std::vector<ulint> sa(n);
for (size_t i = 0; i < n; ++i) {
    sa[n - i - 1] = phi.SA();
    pos = pos.Phi(pos);
}
```

### LF Permutation using RunPerm without RunPermLF
```cpp
// TEXT: GATTACATGATTACATAGATTACATT$
// BWT:  TTTTTCCCGGGAAAT$ATTTTAAAAAA
// RLBWT: TCGAT$ATA
std::vector<std::array<ulint, 1>> bwt_heads = {{'T'},{'C'},{'G'},{'A'},{'T'},{ 0 },{'A'},{'T'},{'A'}};
std::vector<ulint> bwt_run_lengths =          {  5  ,  3  ,  3  ,  3  ,  1  ,  1  ,  1  ,  4  ,  6  };
std::vector<ulint> lf_permutations =          { 17  , 11  , 14  ,  1  , 22  ,  0  ,  4  , 23  ,  5  };
size_t domain = 27;

enum class RunData {
  BWT_CHAR,
  COUNT
};
RunPerm<RunData> rp(bwt_run_lengths, lf_permutations, domain, bwt_heads);
```

## Citation
If using RunPerm or describing length capping in an academic context, please cite:
>Brown, N. K., & Langmead, B. (2026). Bounding the Average Move Structure Query for Faster and Smaller RLBWT Permutations. arXiv. [doi:10.1101/2024.10.29.620953](https://arxiv.org/abs/2602.11029)
