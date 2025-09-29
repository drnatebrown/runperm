#ifndef _PACKED_VECTOR_ALIGNED_HPP
#define _PACKED_VECTOR_ALIGNED_HPP

#include "common.hpp"
#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <iostream>

// Column values are rounded to the nearest byte
template <typename Columns>
class PackedVectorAligned {
public:
    using word_t = uchar;

    // We read ulint at a time, this ensures we never need to read more than one ulint
    constexpr static uchar max_width = NUM_BITS(ulint); // should be 64 bits
    constexpr static size_t NUM_COLS = static_cast<size_t>(Columns::NUM_COLS);

    PackedVectorAligned() = default;
    PackedVectorAligned(const ulint num_rows, const std::array<uchar, NUM_COLS>& widths) {
        PackedVectorAligned::num_rows = num_rows;
        PackedVectorAligned::widths = widths;

        init();
    }

    PackedVectorAligned(PackedVectorAligned&& other) noexcept = default;
    PackedVectorAligned& operator=(PackedVectorAligned&& other) noexcept = default;
    PackedVectorAligned(const PackedVectorAligned& other) = default;
    PackedVectorAligned& operator=(const PackedVectorAligned& other) = default;
    ~PackedVectorAligned() = default;

    template <Columns Col>
    ulint get(size_t i) const {
        assert(i < num_rows);

        size_t byte_pos = get_row_start_byte(i) + byte_offsets[static_cast<size_t>(Col)];
        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        return extract_bits(bits, masks_extract[static_cast<size_t>(Col)]);
    }      

    template <Columns Col>
    void set(size_t i, ulint val) {
        assert(i < num_rows);
        assert(val < POW2(BYTES_TO_BITS(byte_widths[static_cast<size_t>(Col)]))); // value must fit in the column

        size_t byte_pos = get_row_start_byte(i) + byte_offsets[static_cast<size_t>(Col)];

        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        write_bits(bits, masks_write[static_cast<size_t>(Col)], val);
        std::memcpy(&data[byte_pos], &bits, sizeof(ulint));
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
        return byte_row_width*num_rows + sizeof(ulint)/sizeof(word_t);
    }

    size_t get_num_rows() const { return num_rows; }
    size_t get_num_cols() const { return NUM_COLS; }
    std::array<uchar, NUM_COLS> get_widths() const { return widths; }
    std::array<uchar, NUM_COLS> get_byte_widths() const { return byte_widths; }

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
    size_t byte_vector_width; // byte width of stored data (actual data size might be larger due to padding)
    size_t byte_row_width; // width of each row in bytes

    std::array<uchar, NUM_COLS> widths; // Intended bit width of each column
    std::array<uchar, NUM_COLS> byte_widths; // Actual byte width of each column
    std::array<uint16_t, NUM_COLS> byte_offsets; // Offset of the first byte of each column in the vector
    std::array<ulint, NUM_COLS> masks_extract; // Mask of the bits of each column for get
    std::array<ulint, NUM_COLS> masks_write; // Mask of the bits of each column for set

    std::vector<word_t> data;

    inline size_t get_row_start_byte(size_t i) const {
        return i*byte_row_width;
    }

    void init() {
        size_t byte_pos = 0;
        for (size_t i = 0; i < NUM_COLS; i++) {
            assert(widths[i] <= max_width);
            // Round up to the nearest byte
            if (widths[i] % NUM_BITS(word_t) != 0) {
                byte_widths[i] = (widths[i] / NUM_BITS(word_t)) + 1;
            } else {
                byte_widths[i] = widths[i] / NUM_BITS(word_t);
            }

            masks_extract[i] = MASK(widths[i]);
            masks_write[i] = ~(masks_extract[i]);
            
            byte_offsets[i] = byte_pos;
            byte_pos += byte_widths[i];
        }
        byte_row_width = byte_pos;
        byte_vector_width = num_rows * byte_row_width;
        data.resize(data_size());
    }

    inline ulint extract_bits(ulint bits, ulint mask) const {
        return bits & mask;
    }

    inline void write_bits(ulint& bits, ulint mask, ulint val) {
        // clear old bits
        bits &= mask;
        // set new bits
        bits |= val;
    }
};

#endif // end of include guard: _PACKED_VECTOR_ALIGNED_HPP