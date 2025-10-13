CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3

HEADERS = $(shell find include -name "*.hpp")

# all: build invert move_build
all: move_test rlbwt_test runperm_test

# Target to build the executable for build.cpp
# build: src/build.cpp $(HEADERS)
# 	$(CXX) $(CXXFLAGS) -o build src/build.cpp

# # Target to build the executable for invert.cpp
# invert: src/invert.cpp $(HEADERS)
# 	$(CXX) $(CXXFLAGS) -o invert src/invert.cpp

# Target to build the test executable for move_build.cpp
move_test: src/move_test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o move_test src/move_test.cpp

rlbwt_test: src/rlbwt_test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o rlbwt_test src/rlbwt_test.cpp

runperm_test: src/runperm_test.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o runperm_test src/runperm_test.cpp

.PHONY: debug
debug: CXXFLAGS += -g -O0
debug: all

# Clean up build files
clean:
	rm -f build invert move_build move_test rlbwt_test runperm_test