# RunPerm

Efficient data structures for representing run-length encoded permutations, featuring move structures and specialized RLBWT implementations for high-performance indexing.

## Overview

RunPerm provides compact representations of run-length encoded permutations using move structures, with specialized implementations for various indexing applications. The library offers three main interfaces:

- **RunPerm**: Core run-length permutation data structure with integrated move structures
- **MovePerm**: Simplified permutation representation using move structures
- **RLBWT**: Run-length Burrows-Wheeler Transform implementations (in progress)

## Features

- **Compact Representation**: Efficient storage of run-length encoded permutations
- **Move Structures**: Advanced data structures for fast navigation and queries
- **Flexible Configuration**: Multiple template parameters for different use cases
- **High Performance**: Optimized for indexing and pattern matching applications
- **Serialization Support**: Built-in serialization/deserialization capabilities

## Quick Start

### Installation

```bash
git clone https://github.com/drnatebrown/runperm.git
cd runperm
make
```

### Basic Usage

#### RunPerm

```cpp
#include "runperm.hpp"

// Create a run-length permutation
std::vector<ulint> lengths = {5, 3, 7};           // Length of each run
std::vector<ulint> permutation = {0, 2, 1};       // Permutation of runs
ulint domain = 15;                                // Total domain size
std::vector<RunData> run_data = {...};           // Run-specific data

// Basic separated representation (default)
RunPermSeperated<MyRunCols> rp(lengths, permutation, domain, run_data);

// Integrated representation for better cache locality
RunPermIntegrated<MyRunCols> rp_integrated(lengths, permutation, domain, run_data);

// With absolute position storage for index lookups
RunPermIntegratedAbsolute<MyRunCols> rp_absolute(lengths, permutation, domain, run_data);
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

#### RLBWT (TODO)

```cpp
#include "rlbwt.hpp"

// Run-length BWT implementations coming soon
// RunPermLF, MoveLF, and other RLBWT specializations
```

## API Reference

### RunPerm

The main class for run-length encoded permutations with move structure integration.

**Template Parameters:**
- `RunColsType`: Type defining run data columns
- `IntegratedMoveStructure`: Whether to integrate run data with move structure (default: false)
- `StoreAbsolutePositions`: Whether to store absolute positions for index lookups (default: false)

**Key Methods:**
- `next()`: Apply permutation
- `next(ulint steps)`: Apply permutation multiple times
- `up()/down()`: Navigate between runs
- `get_position()/set_position()`: Position management
- `get<Col>()`: Access run data columns
- `pred<Col>(val)/succ<Col>(val)`: Search operations

### MovePerm

Simplified permutation representation using move structures.

**Template Parameters:**
- `StoreAbsolutePositions`: Whether to store absolute positions (default: false)

**Key Methods:**
- `next()`: Apply permutation
- `up()/down()`: Navigate between intervals
- `size()`: Get permutation size
- `move_runs()`: Get number of move structure runs

### MoveStructure

Core move structure implementation for efficient navigation.

**Available Types:**
- `MoveStructureTbl`: Table-based implementation
- `MoveStructureVec`: Vector-based implementation
- `MoveStructureTblIdx`/`MoveStructureVecIdx`: Indexed variants

## Configuration Options

### Default Settings

```cpp
// Default configuration constants
constexpr bool DEFAULT_INTEGRATED_MOVE_STRUCTURE = false;
constexpr bool DEFAULT_STORE_ABSOLUTE_POSITIONS = false;
```

### Run Data Columns

Define custom run data types by implementing a `RunCols` structure:

```cpp
struct MyRunCols {
    static constexpr size_t COUNT = 3;
    enum { CHAR_COUNT, FREQ, OFFSET };
};
```

## Building

### Compilation

```bash
# Standard build
make

# Debug build
make clean
make debug

# Run tests
make test
```

### Dependencies

- C++17 or later
- Standard library support for `std::optional`, `std::array`, etc.

## Examples

### Basic Permutation Navigation

```cpp
#include "runperm.hpp"

// Create permutation
RunPermSeperated<MyRunCols> rp(lengths, permutation, domain, run_data);

// Navigate
rp.first();                    // Move to first run
rp.next();                     // Apply permutation once
rp.next(5);                    // Apply permutation 5 times
rp.down();                     // Move to next run

// Query
auto pos = rp.get_position();
ulint size = rp.size();
ulint runs = rp.move_runs();
```

### Search Operations

```cpp
// Find predecessor/successor
auto pred_pos = rp.pred<MyRunCols::CHAR_COUNT>(target_value);
auto succ_pos = rp.succ<MyRunCols::CHAR_COUNT>(target_value);

if (pred_pos) {
    rp.set_position(*pred_pos);
    // Process found position
}
```

### Serialization

```cpp
// Save to file
std::ofstream out("permutation.dat");
rp.serialize(out);

// Load from file
std::ifstream in("permutation.dat");
RunPermSeperated<MyRunCols> loaded_rp;
loaded_rp.load(in);
```

## Performance Considerations

- **Integrated vs Separated**: Integrated move structures offer better cache locality but may use more memory
- **Absolute vs Relative Positions**: Absolute positions enable faster index lookups but increase memory usage
- **Run Data Size**: Larger run data structures may impact cache performance

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

## License

[License information - update as needed]

## References

- Move structures for efficient permutation representation
- Run-length encoding for compact data structures
- Burrows-Wheeler Transform applications in bioinformatics