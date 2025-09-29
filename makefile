CXX = g++
CXXFLAGS = -std=c++17 -Iinclude -O3

HEADERS = $(wildcard include/*.hpp) $(wildcard include/ds/*.hpp)

# all: build invert move_build
all: move_build

# Target to build the executable for build.cpp
# build: src/build.cpp $(HEADERS)
# 	$(CXX) $(CXXFLAGS) -o build src/build.cpp

# # Target to build the executable for invert.cpp
# invert: src/invert.cpp $(HEADERS)
# 	$(CXX) $(CXXFLAGS) -o invert src/invert.cpp

# Target to build the test executable for move_build.cpp
move_build: src/move_build.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o move_build src/move_build.cpp

.PHONY: debug
debug: CXXFLAGS += -g -O0
debug: all

# Clean up build files
clean:
	rm -f build invert move_build