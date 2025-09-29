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

template <typename Columns>
class PackedVector {
public:
    using word_t = uchar;

    // We read ulint at a time, this ensures we never need to read more than one ulint
    // should be 57 bits for 64 bit ulint and 8 bit word_t
    constexpr static uchar max_width = NUM_BITS(ulint) - (NUM_BITS(word_t) - 1);
    constexpr static size_t NUM_COLS = static_cast<size_t>(Columns::NUM_COLS);

    PackedVector() = default;
    PackedVector(const ulint num_rows, const std::array<uchar, NUM_COLS>& widths) {
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
        return extract_bits(bits, pos.offset, masks_extract[static_cast<size_t>(Col)]);
    }      

    template <Columns Col>
    void set(size_t i, ulint val) {
        assert(i < num_rows);
        assert(val < POW2(widths[static_cast<size_t>(Col)])); // value must fit in the column

        bit_pos pos(get_row_start(i) + offsets[static_cast<size_t>(Col)]);

        ulint bits = 0;
        std::memcpy(&bits, &data[pos.chunk], sizeof(ulint));
        write_bits(bits, pos.offset, masks_write[pos.offset][static_cast<size_t>(Col)], val);
        std::memcpy(&data[pos.chunk], &bits, sizeof(ulint));
    }

    template<size_t... Indices>
    void set_row(size_t i, const std::array<ulint, NUM_COLS>& values, std::index_sequence<Indices...>) {
        (set<static_cast<Columns>(Indices)>(i, values[Indices]), ...);
    }
    void set_row(size_t i, const std::array<ulint, NUM_COLS>& values) {
        set_row(i, values, std::make_index_sequence<NUM_COLS>{});
    }

    template<size_t... Indices>
    std::array<ulint, NUM_COLS> get_row(size_t i, std::index_sequence<Indices...>) const { 
        return {get<static_cast<Columns>(Indices)>(i)...};
    }
    std::array<ulint, NUM_COLS> get_row(size_t i) const { 
        return get_row(i, std::make_index_sequence<NUM_COLS>{});
    }

    size_t size() const { return num_rows; }
    size_t data_size() const { 
        return CEIL_DIV(vector_width, NUM_BITS(word_t)) + sizeof(ulint)/sizeof(word_t);
    }

    size_t get_num_rows() const { return num_rows; }
    size_t get_num_cols() const { return NUM_COLS; }
    std::array<uchar, NUM_COLS> get_widths() const { return widths; }

    // size_t offsets() const { return offsets; }
    // size_t row_width() const { return row_width; }
    // size_t vector_width() const { return vector_width; }

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
        in.read((char *)data.data(), data_size() * sizeof(word_t));
    }

private:
    size_t num_rows;
    size_t vector_width; // bit width of stored data (actual data size might be larger due to padding)
    size_t row_width; // width of each row in bits

    std::array<uchar, NUM_COLS> widths; // Bit width of each column
    std::array<uint16_t, NUM_COLS> offsets; // Offset of the first bit of each column in the vector
    std::array<ulint, NUM_COLS> masks_extract; // Mask of the bits of each column for get
    std::array<std::array<ulint, NUM_COLS>, NUM_BITS(word_t)> masks_write; // The mask for each offset to reset the bits

    std::vector<word_t> data;

    void init() {
        size_t bit_pos = 0;
        for (size_t i = 0; i < NUM_COLS; i++) {
            assert(widths[i] <= max_width);
            offsets[i] = bit_pos;
            masks_extract[i] = MASK(widths[i]);
            for (size_t j = 0; j < NUM_BITS(word_t); j++) {
                masks_write[j][i] = ~(masks_extract[i] << j);
            }
            bit_pos += widths[i];
        }
        row_width = bit_pos;
        vector_width = num_rows * row_width;
        data.resize(data_size());
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

    inline ulint extract_bits(ulint bits, uchar start, ulint mask) const {
        return (bits >> start) & mask;
    }

    inline void write_bits(ulint& bits, uchar start, ulint mask, ulint val) {
        // clear old bits
        bits &= mask;
        // set new bits
        bits |= (val << start);
    }

};

#endif // end of include guard: _PACKED_VECTOR_HPP