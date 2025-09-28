#ifndef _PACKED_VECTOR_HPP
#define _PACKED_VECTOR_HPP

#include "common.hpp"
#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <iostream>

template <typename Columns, typename word_t = uchar>
class PackedVector {
public:
    // We read ulint at a time, this ensures we never need to read more than one ulint
    constexpr static uchar max_width = NUM_BITS(ulint) - 7; // should be 57 bits
    constexpr static size_t NumCols = static_cast<size_t>(Columns::NUM_COLS);

    PackedVector() = default;
    PackedVector(const ulint num_rows, const std::array<uchar, NumCols>& widths) {
        PackedVector::num_rows = num_rows;
        PackedVector::widths = widths;

        init();
    }

    PackedVector(PackedVector&& other) noexcept = default;
    PackedVector& operator=(PackedVector&& other) noexcept = default;
    PackedVector(const PackedVector& other) = default;
    PackedVector& operator=(const PackedVector& other) = default;
    ~PackedVector() = default;

    template <Columns Col>
    ulint get(size_t i) const {
        assert(i < num_rows);

        bit_pos pos(get_row_start(i) + offsets[static_cast<size_t>(Col)]);
        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        return extract_bits(bits, pos.offset, widths[static_cast<size_t>(Col)]);
    }      

    template <Columns Col>
    void set(size_t i, ulint val) {
        assert(i < num_rows);
        assert(val < (1ULL << widths[static_cast<size_t>(Col)])); // value must fit in the column

        bit_pos pos(get_row_start(i) + offsets[static_cast<size_t>(Col)]);

        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        write_bits(bits, pos.offset, widths[static_cast<size_t>(Col)], val);
        std::memcpy(&data[pos.chunk], &bits, sizeof(ulint));
    }

    template<size_t... Indices>
    void set_row(size_t i, const std::array<ulint, NumCols>& values, std::index_sequence<Indices...>) {
        (set<static_cast<Columns>(Indices)>(i, values[Indices]), ...);
    }
    void set_row(size_t i, const std::array<ulint, NumCols>& values) {
        set_row(i, values, std::make_index_sequence<NumCols>{});
    }

    template<size_t... Indices>
    std::array<ulint, NumCols> get_row(size_t i, std::index_sequence<Indices...>) const { 
        return {get<static_cast<Columns>(Indices)>(i)...};
    }
    std::array<ulint, NumCols> get_row(size_t i) const { 
        return get_row(i, std::make_index_sequence<NumCols>{});
    }

    size_t size() const { return num_rows; }
    size_t get_data_size() const { return data_size; }
    size_t get_num_rows() const { return num_rows; }
    size_t get_num_cols() const { return NumCols; }
    size_t get_row_width() const { return row_width; }

    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        out.write((char *)&num_rows, sizeof(num_rows));
        written_bytes += sizeof(num_rows);

        out.write((char *)widths.data(), sizeof(widths));
        written_bytes += sizeof(widths);

        
        const char* data_mem = reinterpret_cast<const char*>(data.data());
        size_t size = data.size() * sizeof(word_t);
        out.write(data_mem, size);
        written_bytes += size;

        return written_bytes;
    }

    void load(std::istream &in) {
        in.read((char *)&num_rows, sizeof(num_rows));
        in.read((char *)widths.data(), sizeof(widths));
        init();
        in.read((char *)data.data(), data_size * sizeof(word_t));
    }

private:
    size_t num_rows;
    size_t vector_width; // bit width of stored data (actual data size might be larger due to padding)
    size_t data_size; // number of words in data vector
    size_t row_width; // width of each row in bits

    std::array<uchar, NumCols> widths; // Bit width of each column
    std::array<uint16_t, NumCols> offsets; // Offset of the first bit of each column in the vector

    std::vector<word_t> data;

    void init() {
        size_t bit_pos = 0;
        for (size_t i = 0; i < NumCols; i++) {
            assert(widths[i] <= max_width);
            offsets[i] = bit_pos;
            bit_pos += widths[i];
        }
        row_width = bit_pos;
        vector_width = num_rows * row_width;
        data_size = CEIL_DIV(vector_width, NUM_BITS(word_t)) + sizeof(ulint)/sizeof(word_t);
        data.resize(data_size);
    }

    struct bit_pos {
        ulint chunk;
        ulint offset;

        bit_pos(size_t bit) {
            chunk = bit / NUM_BITS(word_t);
            offset = bit % NUM_BITS(word_t);
        }
    };

    inline size_t get_row_start(size_t i) const {
        return i*row_width;
    }

    inline ulint extract_bits(ulint bits, uchar start, uchar width) const {
        return (bits >> start) & ((1ULL << width) - 1);
    }

    inline void write_bits(ulint& bits, uchar start, uchar width, ulint val) {
        // clear old bits
        bits &= ~(((1ULL << width) - 1) << start);
        // set new bits
        bits |= (val << start);
    }

};

#endif // end of include guard: _PACKED_VECTOR_HPP