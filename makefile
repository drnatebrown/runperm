VERSION = 0.1.3

CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3

HEADERS = $(shell find include -name "*.hpp")

BUILD_DIR = build
UNIT_BUILD_DIR = $(BUILD_DIR)/unit
INTEGRATION_BUILD_DIR = $(BUILD_DIR)/integration
BENCH_BUILD_DIR = $(BUILD_DIR)/bench

# High-level targets
# - all: build all unit/integration tests + benchmarks + examples
# - test: build unit/integration tests and run them

UNIT_TESTS = $(UNIT_BUILD_DIR)/packed_vector_test \
             $(UNIT_BUILD_DIR)/alphabet_test \
             $(UNIT_BUILD_DIR)/move_table_test \
             $(UNIT_BUILD_DIR)/move_splitting_test \
             $(UNIT_BUILD_DIR)/move_structure_test \
             $(UNIT_BUILD_DIR)/moveperm_test \
             $(UNIT_BUILD_DIR)/runperm_test \
             $(UNIT_BUILD_DIR)/rlbwt_row_test \
             $(UNIT_BUILD_DIR)/rlbwt_structure_test
INTEGRATION_TESTS = $(INTEGRATION_BUILD_DIR)/rlbwt_test \
                    $(INTEGRATION_BUILD_DIR)/move_structure_test \
                    $(INTEGRATION_BUILD_DIR)/moveperm_test \
                    $(INTEGRATION_BUILD_DIR)/runperm_test
BENCH_TESTS = $(BENCH_BUILD_DIR)/move_bench \
              $(BENCH_BUILD_DIR)/runperm_bench \
              $(BENCH_BUILD_DIR)/rlbwt_bench

all: $(UNIT_TESTS) $(INTEGRATION_TESTS) $(BENCH_TESTS) examples

test: $(UNIT_TESTS) $(INTEGRATION_TESTS)
	$(UNIT_BUILD_DIR)/packed_vector_test
	$(UNIT_BUILD_DIR)/alphabet_test
	$(UNIT_BUILD_DIR)/move_table_test
	$(UNIT_BUILD_DIR)/move_splitting_test
	$(UNIT_BUILD_DIR)/move_structure_test
	$(UNIT_BUILD_DIR)/moveperm_test
	$(UNIT_BUILD_DIR)/runperm_test
	$(UNIT_BUILD_DIR)/rlbwt_row_test
	$(UNIT_BUILD_DIR)/rlbwt_structure_test
	$(INTEGRATION_BUILD_DIR)/rlbwt_test
	$(INTEGRATION_BUILD_DIR)/move_structure_test
	$(INTEGRATION_BUILD_DIR)/moveperm_test
	$(INTEGRATION_BUILD_DIR)/runperm_test
	@echo "================================================="
	@echo "     All unit and integration tests passed"
	@echo "================================================="

examples: examples.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Unit tests (header-only data structures)
$(UNIT_BUILD_DIR)/packed_vector_test: ./tests/unit/ds/packed_vector_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/alphabet_test: ./tests/unit/ds/alphabet_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_table_test: ./tests/unit/move/move_table_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_splitting_test: ./tests/unit/move/move_splitting_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/move_structure_test: ./tests/unit/move/move_structure_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/moveperm_test: ./tests/unit/runperm/moveperm_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/runperm_test: ./tests/unit/runperm/runperm_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/rlbwt_row_test: ./tests/unit/rlbwt/rlbwt_row_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(UNIT_BUILD_DIR)/rlbwt_structure_test: ./tests/unit/rlbwt/rlbwt_structure_test.cpp $(HEADERS)
	mkdir -p $(UNIT_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Integration-style tests that exercise larger rlbwt/runperm flows
$(INTEGRATION_BUILD_DIR)/rlbwt_test: ./tests/integration/rlbwt_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/move_structure_test: ./tests/integration/move_structure_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/moveperm_test: ./tests/integration/moveperm_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(INTEGRATION_BUILD_DIR)/runperm_test: ./tests/integration/runperm_test.cpp $(HEADERS)
	mkdir -p $(INTEGRATION_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

# Benchmarks (not run by default in `make test`)
$(BENCH_BUILD_DIR)/move_bench: ./tests/benchmark/move_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/runperm_bench: ./tests/benchmark/runperm_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(BENCH_BUILD_DIR)/rlbwt_bench: ./tests/benchmark/rlbwt_bench.cpp $(HEADERS)
	mkdir -p $(BENCH_BUILD_DIR)
	$(CXX) $(CXXFLAGS) -o $@ $<

.PHONY: debug
debug: CXXFLAGS += -g -O0
debug: all

# Clean up build files
clean:
	rm -rf $(BUILD_DIR)
	rm -f examples