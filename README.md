# RunPerm

Efficient, header-only data structures for representing run-length encoded permutations, featuring move structures and specialized RLBWT implementations for high-performance indexing.

## Overview

RunPerm provides compact representations of run-length encoded permutations using move structures, with specialized implementations for various indexing applications. The library offers three main interfaces:

- **RunPerm**: Core run-length permutation data structure with integrated move structures
- **MovePerm**: Simplified permutation representation using move structures
- **RLBWT**: Run-length Burrows-Wheeler Transform implementations (LF/FL and Phi/InvPhi) built on RunPerm and MovePerm

## Features

- **Compact Representation**: Efficient storage of run-length encoded permutations
- **Move Structures**: Advanced data structures for fast navigation and queries
- **RLBWT Specializations**:
  - LF/FL navigation with character access
  - Phi/InvPhi navigation with SA sampling
  - Helpers to derive Phi/InvPhi structures from run-length BWT
- **Flexible Configuration**: Multiple template parameters for different use cases
- **High Performance**: Optimized for indexing and pattern matching applications
- **Serialization Support**: Built-in serialization/deserialization capabilities

## Quick Start

### Installation

This library is header-only. There is no library artifact to build or install. Just add `include/` to your compiler's include paths and include `runperm.hpp`, `move.hpp`, or `rlbwt.hpp` as needed.

The provided `makefile` only builds example/test executables:
- `move_test`
- `rlbwt_test`
- `runperm_test`

### Basic Usage

#### RunPerm

```cpp
#include "runperm.hpp"

// Create a run-length permutation
std::vector<ulint> lengths = {5, 3, 7};           // Length of each run
std::vector<ulint> permutation = {0, 2, 1};       // Permutation of runs
ulint domain = 15;                                // Total domain size

// Your run data (example)
struct MyRunCols {
    static constexpr size_t COUNT = 3;
    enum { CHAR_COUNT, FREQ, OFFSET };
};
using RunData = std::array<ulint, MyRunCols::COUNT>;
std::vector<RunData> run_data(/* lengths.size() */);

// Basic construction
RunPerm<MyRunCols> rp(lengths, permutation, domain, run_data);

```

#### MovePerm

```cpp
#include "runperm.hpp"

// Create from permutation vector
std::vector<ulint> permutation = {2, 0, 1, 3};
MovePermRelative mp(permutation);

// Or from run structure
MovePermAbsolute mp_abs(lengths, permutation, domain);
```

#### RLBWT: LF/FL
```cpp
#include "rlbwt.hpp"

// BWT = TTTTTCCCGGGAAAT$ATTTTAAAAAA
// RLBWT as run heads (characters) and run lengths, with 1 as terminal
std::vector<uchar> bwt_heads = {'T','C','G','A','T', 1 ,'A','T','A'};
std::vector<ulint> bwt_run_lengths = {5,3,3,3,1,1,1,4,6};

// LF using Move-only interface
MoveLF<> move_lf(bwt_heads, bwt_run_lengths);

// LF with RunPerm + run data columns
enum class RunCols { VAL_1, VAL_2, COUNT };
using LFRunData = std::array<ulint, (size_t)RunCols::COUNT>;
std::vector<LFRunData> lf_run_data(bwt_heads.size()); // Insert with some type of data
RunPermLF<RunCols> runperm_lf(bwt_heads, bwt_run_lengths, lf_run_data);

// FL (symmetrically)
MoveFL<> move_fl(bwt_heads, bwt_run_lengths);
...
```

#### RLBWT: Phi/InvPhi

```cpp
#include "rlbwt.hpp"

// Build Phi structure from RLBWT
auto [phi_lengths, phi_interval_perm, domain] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);

// RunPerm Phi (always absolute positions)
enum class PhiCols { VAL_1, VAL_2, COUNT };
using PhiRunData = std::array<ulint, (size_t)PhiCols::COUNT>;
std::vector<PhiRunData> phi_run_data(phi_lengths.size());
RunPermPhi<PhiCols> phi(phi_lengths, phi_interval_perm, domain, phi_run_data);

// Move-only Phi
auto [phi_lens2, phi_perm2, domain2] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
MovePhi mphi(phi_lens2, phi_perm2, domain2);

// Build InvPhi similarly...
auto [invphi_lengths, invphi_interval_perm, domain_inv] = rlbwt_to_invphi(bwt_heads, bwt_run_lengths);
RunPermInvPhi<PhiCols> invphi(invphi_lengths, invphi_interval_perm, domain_inv, phi_run_data);
...

```

## Interface

### RunPerm

The main class for run-length encoded permutations with move structure integration.

- **Template Parameters**:
  - `RunColsType`: Type defining run data columns
  - `IntegratedMoveStructure`: integrate run data with move structure (default: false)
  - `StoreAbsolutePositions`: store absolute positions for index lookups (default: false)
- **Key Methods**:
  - `next()`, `next(ulint steps)`
  - `up()`, `down()`, `first()`, `last()`
  - `get_position()`, `set_position(pos)`
  - `get<Col>()`, `get<Col>(i)`, `get_length()`, `get_length(i)`
  - `pred<Col>(val)`, `succ<Col>(val)`
  - `size()`, `move_runs()`, `permutation_runs()`
  - `serialize(os)`, `load(is)`

### MovePerm

Simplified permutation representation using move structures.

- **Template Parameters**:
  - `StoreAbsolutePositions` (default: false)
- **Key Methods**:
  - `next()`, `up()`, `down()`, `first()`, `last()`
  - `size()`, `move_runs()`, `permutation_runs()`
  - `get_length()`, `get_length(i)`
  - `serialize(os)`, `load(is)`

### RLBWT: LF/FL

- **Types**:
  - `RunPermLF<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType>`
  - `MoveLF<StoreAbsolutePositions, AlphabetType>`
  - `RunPermFL<RunColsType, IntegratedMoveStructure, StoreAbsolutePositions, AlphabetType>`
  - `MoveFL<StoreAbsolutePositions, AlphabetType>`
  - Convenience aliases: `RunPermLFSeperated`, `RunPermLFIntegrated`, `RunPermLFSeperatedAbsolute`, `RunPermLFIntegratedAbsolute` and similarly for FL.
- **Construction**:
  - From `std::vector<uchar> bwt_heads` and `std::vector<ulint> bwt_run_lengths`
- **Key Methods**:
  - `LF()`, `LF(ulint steps)`; `FL()`, `FL(ulint steps)`
  - `get_character()` (at current position)
  - `get_character(run_index)` (character for a given run)
  - All standard RunPerm/MovePerm navigation and queries
- **Alphabet**:
  - Default `AlphabetType = Nucleotide`. Internally maps inputs and handles special `TERMINATOR` and `SEPARATOR`.

### RLBWT: Phi/InvPhi and Helpers

- **Helpers**:
  - `rlbwt_to_phi(bwt_heads, bwt_run_lengths)` -> `(lengths, interval_permutation, domain)`
  - `rlbwt_to_phi(bwt_heads, bwt_run_lengths, lf)` overload to reuse an LF instance
  - `rlbwt_to_invphi(bwt_heads, bwt_run_lengths)` and LF overload similarly
- **Types**:
  - `RunPermPhi<RunColsType, IntegratedMoveStructure, TableType>`
  - `MovePhi` (absolute positions)
  - `RunPermInvPhi<RunColsType, IntegratedMoveStructure, TableType>`
  - `MoveInvPhi` (absolute positions)
- **Key Methods**:
  - `Phi()`, `Phi(ulint steps)`, `SA()` for Phi
  - `InvPhi()`, `InvPhi(ulint steps)`, `SA()` for InvPhi
  - Always uses absolute positions to expose `idx` via `SA()`

### Dependencies

- C++17 or later
- No external dependencies; header-only usage.

## Examples

### Basic Permutation Navigation

```cpp
#include "runperm.hpp"

RunPerm<MyRunCols> rp(lengths, permutation, domain, run_data);

// Navigate
rp.first();
rp.next();
rp.next(5);
rp.down();

// Query
auto pos = rp.get_position();
ulint size = rp.size();
ulint runs = rp.move_runs();
```

### RLBWT LF Text Recovery

```cpp
#include "rlbwt.hpp"

MoveLF<> lf(bwt_heads, bwt_run_lengths);
lf.first();
std::string recovered(text.size(), '\0');
for (size_t i = 1; i < lf.size(); ++i) {
    recovered[text.size() - i] = (char)lf.get_character();
    lf.LF();
}
```

### Phi SA Recovery

```cpp
#include "rlbwt.hpp"

auto [phi_lens, phi_perm, n] = rlbwt_to_phi(bwt_heads, bwt_run_lengths);
MovePhi phi(phi_lens, phi_perm, n);
phi.last();
phi.Phi();
std::vector<ulint> sa(n);
for (size_t i = 0; i < n; ++i) {
    sa[n - i - 1] = phi.SA();
    phi.Phi();
}
```

## Performance Considerations

- **Integrated vs Separated**: Integrated move structures offer better cache locality but may use more memory
- **Absolute vs Relative Positions**: Absolute positions enable faster index lookups but increase memory usage
- **Run Data Size**: Larger run data structures may impact cache performance
- **RLBWT Mapping**: Alphabet mapping and special symbols handling can affect constant factors
