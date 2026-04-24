# Orbit
![GitHub release](https://img.shields.io/github/v/release/drnatebrown/orbit)

$O(r)$-space bitpacked (Orbit) move structures for runny permutations! 

Flexible plug-and-play header library implementing compact data structures for representing run-length encoded permutations. Since permutations decompose into cycles, navigation becomes a literal orbit, stepping through positions as if circling an encapsulated world. Engineered for compactness, Orbit features $O(r)$-space bitpacked move structures and specialized RLBWT implementations for high-performance indexing. Can be used as a foundation for further applications, providing storage and retrieval of user fields alongside the permutation intervals. Further, powerful splitting optimizations such as length capping and balancing are implemented. 

Described in greater detail at [https://arxiv.org/abs/2602.11029](https://arxiv.org/abs/2602.11029).

## Quick Start

`orbit::permutation` is the generic implementation allowing storage of user defined fields alongside a runny permutation.

Consider this runny permutation example:

<img width="598" height="572" alt="Image" src="https://github.com/user-attachments/assets/18933a68-b858-40eb-bfb3-f5fd758557d6" />

We pass the lengths of contiguously permuted intervals ($S_\ell$) and the permutation of their first values (images, $S_\pi$):
```cpp
// Create a run-length permutation
// Length of contiguously permuted intervals
std::vector<unsigned long int> lengths = {2, 3, 1,  2, 2,  1, 1,  1, 3};
// Permutation at head of each interval
std::vector<unsigned long int> images  = {1, 9, 3, 12, 4, 14, 0, 15, 6};   
```

These can be used to build a permutation object, which permits efficient navigation:
```cpp
#include "orbit/permutation.hpp"

orbit::permutation perm(lengths, images);
auto pos = perm.first();
pos = perm.next(pos);
...
```

## Overview

Orbit provides compact representations of run-length encoded permutations, with specialized implementations for various indexing applications. The core permutation logic is stored in a move structure. The library offers four main components:

- **Permutation**: Core run-length permutation data structure allowing user defined fields stored alongside permutation intervals.
- **RLBWT**: Specialized RLBWT permutations (LF/FL and $\phi$ / $\phi^{-1}$) built on permutation features above.
- **Move Structure**: Implementation of the "move data structure", foundation of above methods.
- **Interval Encoding**: Decomposition into contigiously permuted intervals, used to construct above methods. Allows independent use of splitting algorithms such as length capping and balancing.

## Features

- **Move Structures**: Advanced data structures for fast and cache efficient navigation of runny permutations.
- **Compact Representation**: Bitpacked, using the minimum fixed width per move structure column and user defined data columns.
- **Advanced Storage** Allows efficient storage and retrieval of additional information which can be bitpacked alongside the permutation for cache efficiency.
- **Splitting Schemes**: Optimization schemes for length capping and balancing providing faster and smaller move structures, running in optimal $O(r)$-time and space.
- **RLBWT Specializations**:
  - LF/FL navigation with character access, constructed in $O(r)$-time and space from RLBWT.
  - $\phi$ / $\phi^{-1}$ navigation with SA retrieval.
  - Helpers to derive $\phi$ / $\phi^{-1}$ permutations from RLBWT in $O(n)$-time and $O(r)$-space.
- **Flexible Configuration**: Multiple template parameters for various move structure representations.

## Guide

### Installation

This library is header-only. Just add `include/` to your compiler's include paths and include `orbit/permutation.hpp`, `orbit/move_structure.hpp`, `orbit/rlbwt.hpp`, or `orbit/interval_encoding.hpp` as needed:

- `permutation.hpp` loads ``orbit::permutation`` and its dependencies.
- `move_structure.hpp` loads the underlying ``orbit::move_structure`` implementation only.
- `rlbwt.hpp` loads specialized classes to represent LF/FL and $\phi$ / $\phi^{-1}$. These are ``orbit::rlbwt::lf_permutation``, ``orbit::rlbwt::fl_permutation``, ``orbit::rlbwt::phi_permutation``, and ``orbit::rlbwt::phi_inv_permutation``.
- ``interval_encoding.hpp`` loads ``orbit:interval_encoding`` which converts a permutation into its run/interval length encoding with option to apply splitting algorithms.

The provided `makefile` only builds example/test executables:

- `make example` builds the `example` executable with helpful use cases.
- `make test` builds and runs all test modules.
- `make bench` builds benchmarking executables.

### Basic Usage

#### Permutation

See quick start above! Included with `orbit/permutation.hpp` where further documentation is found. Some small examples are given in the provided `examples/example.cpp` file. Although the example passed the run-length encoded permutation, the `orbit::permutation` class can also accept the raw permutation $\pi[0..n-1]$ itself.

#### Data Columns

We can also define a number of columns of additional data to store alongside the permutation intervals. This is made easier with a macro:

```cpp
DEFINE_ORBIT_COLUMNS(data_cols, VAL1, VAL2)
// The DEFINE_COLUMNS(enum_name, ...) macro above is equivalent to:
enum class data_cols {
     VAL1,
     VAL2,
     COUNT
};
// The COUNT enumerator is automatically added by the DEFINE_COLUMNS macro,
// but must be included manually in the explicit enum definition.
```

The actual data is then stored in a vector of tuples. The length of this vector must match the number of permutation intervals.

```cpp
// Defines an array of unsigned long integers of length COUNT
using data_tuple = orbit::columns_tuple<data_cols>;
// This is equivalent to
using data_tuple = std::array<unsigned long int, static_cast<size_t>(data_cols::COUNT>)

// We must have some data values for each row
std::vector<data_tuple> run_data(lengths.size());
```

Given this information, we pass the user defined run columns to construct a runny permutation:

```cpp
orbit::permutation<data_cols> meta_perm(lengths, images, run_data);
```

Finally, we use the position object to navigate the permutation and access data:
```cpp
using position = typename orbit::permutation<data_cols>::position
position pos = meta_perm.first(); // start from 0
unsigned long int some_data = meta_perm.get<data_cols::VAL1>(pos);
pos = meta_perm.next(pos); // move one permutation step
unsigned long int other_data = meta_perm.get<data_cols::VAL2>(pos);
```

#### RLBWT: LF/FL

By including `orbit/rlbwt.hpp` users can use specialized methods designed for permutations based on the RLBWT such as LF/FL for pattern matching and text extraction. The library can also first parse a BWT into its RLBWT.

```cpp
#include "orbit/rlbwt.hpp"

// BWT = TTTTTCCCGGGAAAT$ATTTTAAAAAA
// RLBWT as run heads (characters) and run lengths, with 0 as terminator
std::vector<uchar> bwt_heads                   = {'T','C','G','A','T', 0 ,'A','T','A'};
std::vector<unsigned long int> bwt_run_lengths = { 5 , 3 , 3 , 3 , 1 , 1 , 1 , 4 , 6 };

orbit::rlbwt::lf_permutation<> lf(bwt_heads, bwt_run_lengths);

// Can also add data columns
enum class data_cols { VAL1, VAL2, COUNT };
...
orbit::rlbwt::lf_permutation<data_cols> meta_lf(bwt_heads, bwt_run_lengths, data);

// FL (symmetrically)
orbit::rlbwt:fl_permutation fl(bwt_heads, bwt_run_lengths);
...
```

#### RLBWT: $\phi$ / $\phi^{-1}$

By including `orbit/rlbwt.hpp`, we can also build the permutations from an RLBWT for $\phi$ and $\phi^{-1}$ needed for locate queries.

```cpp
#include "orbit/rlbwt.hpp"

orbit::rlbwt::phi_permutation<> phi(bwt_heads, bwt_run_lengths);

// Build phi_inv similarly... Can use user data here too!
DEFINE_ORBIT_COLUMNS(phi_data_cols, VAL1, VAL2);
...
orbit::rlbwt::phi_inv_permutation<phi_data_cols> phi_inv(bwt_heads, bwt_run_lengths, phi_data);
...

```

## Interface

The main class interfaces are summarized by the following components, where `ulint` refers to `unsigned long int`:

- **Constructors**
  - No User Data:
    ```cpp
    orbit::permutation<>(std::vector<ulint> permutation, split_params sp = split_params())
    orbit::permutation<>(std::vector<ulint> lengths, std::vector<ulint> images, split_params sp = split_params())
    orbit::permutation<>(interval_encoding enc)
    ```
  - With User Data (see advanced usage for effect of split_params)
    ```cpp
    using tuple = orbit::columns_tuple<data_t>
    orbit::permutation<data_t>(std::vector<ulint> lengths, std::vector<ulint> images, std::vector<tuple> data)
    orbit::permutation<data_t>(std::vector<ulint> lengths, std::vector<ulint> images, split_params sp, std::vector<tuple> data)
    orbit::permutation<data_t>(interval_encoding enc, std::vector<tuple> data)
    ```
  - `orbit::interval_encoding` offers many alternative input parameters apart from lengths/images.
- **Template Parameters**:
  - `data_columns_t`: Type defining user run data columns (default: empty_data_columns). Alias move_permutation implicitly sets this to empty_data_columns.
  - `integrated_move_structure`: integrate run data bitpacked alongside move structure if true, or stored bitpacked in its own table if false (default: false)
  - `store_absolute_positions`: store absolute positions for index lookups, rather than just interval/offset paits (default: false)
- **Key Methods**:
  - `next(pos)`, `next(pos, ulint steps)`
  - `up(pos)`, `down(pos)`, `first()`, `last()`
  - `get<col>(pos)`, `get<col>(i)`, `get_length(pos)`, `get_length(i)`
  - `pred<col>(pos, val)`, `succ<col>(pos, val)`
  - `domain()`, `runs()`, `intervals()`
  - `serialize(os)`, `load(is)`
    
The public API simplifies template parameters and methods, see the internal implementation for advanced flexibility for building move structure types (discussed below in advanced usage).

## Performance Considerations

- **Integrated vs Separated**: Integrating user data alongside the move structure offer better cache locality but may cause slower move queries since navigating the data structure requires loading larger entries.
- **Absolute vs Relative Positions**: Absolute positions enable full permutation positional information but increases memory usage. The space usage for a runny permutation of $r$ runs over domain $n$ is approximately:
  - **Absolute**: $r \log r + 2 r \log n$ bits
  - **Relative**: $r \log r + 2 r \log \frac{n}{r}$ bits
- **Length Capping**: Can greatly reduce the size of the data structure, especially for relative positions. Also gives amortized guarantees and practical speed up when tuned correctly. Where length capping factor is $c$, splits any runs/intervals longer than $c$ times the average run length of the original permutation. Takes $O(r)$-time and space ([see here](https://arxiv.org/abs/2602.11029)). Where default splitting parameters are used, $c$ is set to 8.
- **Balancing**: Where $\alpha$ is the balancing factor, guarantees less than $2\alpha$ complexity for a single permutation step. However, can increase the size of the data structures due to splitting intervals. Length capping and balancing often work well together. This implements the first $O(r)$-time and space algorithm ([see here](https://arxiv.org/abs/2603.22147)). Where default splitting parameters are used, $\alpha$ is set to 16.

## Advanced Usage

#### Interval Splitting
When we build a permutation object without data columns, default splitting is used. Where a run refers to a maximal contigiously permuted subsequence, an interval is a potentially non-maximal (sub-run) found by splitting runs/intervals. Given this example:

```cpp
std::vector<unsigned long int> lengths = {2, 3, 1,  2, 2,  1, 1,  1, 3};
std::vector<unsigned long int> images  = {1, 9, 3, 12, 4, 14, 0, 15, 6};
orbit::permutation<> perm(lengths, images);
// This is equivalent to
orbit::permutation<> perm(lengths, images, orbit::split_params());
```
The underlying representation may have more intervals than original runs due to under the hood optimizations which make the Orbit library powerful. Thus, we find that these values may be different:
```cpp
perm.runs();
perm.intervals();
```

When a user specifies data columns, splitting optimizations are turned off by default:
```cpp
std::vector<unsigned long int> lengths = {2, 3, 1,  2, 2,  1, 1,  1, 3};
std::vector<unsigned long int> images  = {1, 9, 3, 12, 4, 14, 0, 15, 6};
DEFINE_ORBIT_COLUMNS(data_columns, VAL1, VAL2);
...
orbit::permutation<data_columns> perm(lengths, images, run_data);
//This is equivalent to
orbit::permutation<data_columns> perm(lengths, images, run_data, orbit::NO_SPLITTING);
```

This is because it is not obvious if user specified data can be preserved when a run/interval is split. If a user specifies some splitting parameters, data is copied on split, meaning both new intervals contain the same data values as the former before splitting.
```cpp
// Same as passing orbit::split_params(), the default
orbit::permutation<data_columns> perm(lengths, images, run_data, orbit::split_params(8.0, 16)); // Length cap, balance with param 8/16 resp.
orbit::permutation<data_columns> perm(lengths, images, run_data, orbit::split_params(8.0, std::nullopt)); // Length cap, no balancing
orbit::permutation<data_columns> perm(lengths, images, run_data, orbit::split_params(std::nullopt, 16)); // Balancing, no length capping
```

If a user desires splitting but data cannot be trivially copied, they should first build the interval encoding which contains the final intervals and their lengths. Users can then manually amend their data values before constructing from the encoding:
```cpp
orbit::interval_encoding enc(lengths, images, orbit::split_params);
std::vector<data_tuple> data(enc.intervals());
for(size_t i = 0; i < enc.intervals(); ++i) {
    size_t interval_length = enc.get_length(i);
    data[i] = // some computation
}

orbit::permutation perm(enc, data);
```
#### RLBWT Interface
The RLBWT permutations LF/FL have extra components (we show use for LF below):
- **RLBWT Specific Constructors**
  - No User Data
    ```cpp
    orbit::rlbwt::lf_permutation<>(std::vector<char> bwt, split_params sp = split_params())
    template<typename container_t>
    orbit::rlbwt::lf_permutation<>(container_t rlbwt, container_t rlbwt_lengths, split_params sp = split_params())
    template<typename rlbwt_interval_encoding_t>
    orbit::rlbwt::lf_permutation<>(rlbwt_interval_encoding_t enc)
    ```
  - With User Data (see advanced usage for effect of split_params)
    ```cpp
    using tuple = orbit::columns_tuple<data_t>
    orbit::rlbwt::lf_permutation<>(container_t rlbwt, container_t rlbwt_lengths, std::vector<tuple> data)
    orbit::rlbwt::lf_permutation<>(container_t rlbwt, container_t rlbwt_lengths, split_params sp, std::vector<tuple> data)
    template<typename rlbwt_interval_encoding_t>
    orbit::rlbwt::lf_permutation<>(rlbwt_interval_encoding_t enc, std::vector<tuple> data)
    ```
- **RLBWT Specific Template Parameters**
  - `alphabet`: when using LF/FL permutations, specify `nucleotide` if using DNA alphabets (with 0 reserved for terminators, 1 reserved for separators). Use `alphabet` otherwise.
- **RLBWT Specific Methods**
  - `get_character(pos)`

#### Template Parameters
Many of the library classes expose a simplified interface, including only the template parameters described above. More advanced template parameters can be accessed by loading the specific implementation of a class from the `orbit/include/internal` filetree. For example, `orbit::permutation_impl` exposes parameters that turn on an exponential search for navigation steps instead of the default linear scan (can be preferred when no splitting optimizations are used) or allows modification of the bitpacking scheme and underlying layout of the move structure data. These specific options are likely to be of interest to an academic community, but not for those looking only for a practical library.

## Examples

### Basic Permutation Navigation

```cpp
#include "orbit/permutation.hpp"

orbit::permutation<data_columns> perm(lengths, images);

// Navigate
auto pos = perm.first();
pos = perm.next(pos); // One permutation step
pos = perm.next(pos, 5); // Five permutation steps
pos = perm.down(pos); // Move down to the next interval
pos = perm.up(pos); // Move up to the previous interval

pos.interval; // Current interval containing position
pos.offset; // Offset within that interval
pos.idx; // If using absolute positions (see performance considerations) the actual position in [0, n)

// Largest position less than or equal to this position containing val in column VAL_1
auto pred_pos = perm.pred<data_columns::VAL_1>(pos, val);
auto succ_pos = perm.succ<data_columns::VAL_2>(pos, val); // Similar, smallest larger or equal to
```

### RLBWT LF Text Recovery

```cpp
#include "orbit/rlbwt.hpp"

orbit::rlbwt::lf_permutation lf(bwt_heads, bwt_run_lengths);
auto pos = lf.first();
std::string recovered(lf.domain(), '\0');
for (size_t i = 1; i < lf.domain(); ++i) {
    recovered[lf.domain() - i] = lf.get_character(pos);
    pos = lf.LF(pos); // or, lf.next(pos)
}
```

### $\phi$ SA Recovery

```cpp
#include "orbit/rlbwt.hpp"

orbit::rlbwt::phi_permutation phi(bwt_heads, bwt_lengths);
auto pos = phi.last();
pos = phi.phi(pos); // or, phi.next(pos)
std::vector<ulint> sa(phi.domain());
for (size_t i = 0; i < phi.domain(); ++i) {
    sa[phi.domain() - i - 1] = phi.SA(); // or, pos.idx
    pos = pos.phi(pos);
}
```

### LF Permutation using Permutation without RLBWT Specializations

```cpp
DEFINE_ORBIT_COLUMNS(bwt_data, BWT_CHAR);
using char_tuple = columns_tuple(bwt_data);

// TEXT:  GATTACATGATTACATAGATTACATT$
// BWT:   TTTTTCCCGGGAAAT$ATTTTAAAAAA
// RLBWT: T    C  G  A  T$AT   A
std::vector<char_tuple> bwt_heads           = {{'T'},{'C'},{'G'},{'A'},{'T'},{ 0 },{'A'},{'T'},{'A'}};
std::vector<ulint> bwt_run_lengths          = {  5  ,  3  ,  3  ,  3  ,  1  ,  1  ,  1  ,  4  ,  6  };
std::vector<ulint> lf_images                = { 17  , 11  , 14  ,  1  , 22  ,  0  ,  4  , 23  ,  5  };

orbit::permutation rlbwt(bwt_run_lengths, lf_images, bwt_heads);
```

## Dependencies

- C++17 or later
- No external dependencies; header-only usage.

## Citation

If using Orbit or describing length capping in an academic context, please cite:

> Brown, N. K., & Langmead, B. (2026). Bounding the Average Move Structure Query for Faster and Smaller RLBWT Permutations. arXiv. [doi.org/10.48550/2602.11029](https://arxiv.org/abs/2602.11029)

If you rely on Orbit for, or describe, optimal-time move structure construction or balancing, please cite:

> Brown, N.K., Sanaullah, A., Zhang, S., & Langmead, B. (2026). Optimal-Time Move Structure Construction. arXiv. [doi.org/10.48550/2603.22147](https://arxiv.org/abs/2603.22147)

