#ifndef _PACKED_VECTOR_ALIGNED_HPP
#define _PACKED_VECTOR_ALIGNED_HPP

#include "orbit/common.hpp"
#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <iterator>

namespace orbit {

// Column values are rounded to the nearest byte
// TODO this would be better if definitely fast, using SIMD to unpack entire rows?
template <size_t num_cols>
class packed_matrix_aligned {
public:
    using word_t = uchar;

    // We read ulint at a time, this ensures we never need to read more than one ulint
    constexpr static uchar max_width = num_bits_type(ulint); // should be 64 bits

    packed_matrix_aligned() = default;
    packed_matrix_aligned(const ulint rows, const std::array<uchar, num_cols>& widths) {
        packed_matrix_aligned::num_rows = rows;
        packed_matrix_aligned::widths = widths;

        init();
    }

    packed_matrix_aligned(packed_matrix_aligned&& other) noexcept = default;
    packed_matrix_aligned& operator=(packed_matrix_aligned&& other) noexcept = default;
    packed_matrix_aligned(const packed_matrix_aligned& other) = default;
    packed_matrix_aligned& operator=(const packed_matrix_aligned& other) = default;
    ~packed_matrix_aligned() = default;

    template<size_t col>
    ulint get(size_t row) const {
        static_assert(col < num_cols, "Column out of bounds");
        assert(row < num_rows);

        size_t byte_pos = get_row_start_byte(row) + byte_offsets[col];
        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        return extract_bits(bits, masks_extract[col]);
    }      

    template<size_t col>
    void set(size_t row, ulint val) {
        static_assert(col < num_cols, "Column out of bounds");
        assert(row < num_rows);
        assert(val < pow2(bytes_to_bits(byte_widths[col]))); // value must fit in the column

        size_t byte_pos = get_row_start_byte(row) + byte_offsets[col];

        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        write_bits(bits, masks_write[col], val);
        std::memcpy(&data[byte_pos], &bits, sizeof(ulint));
    }

    template<size_t... col>
    void set_row(size_t row, const std::array<ulint, num_cols>& values, std::index_sequence<col...>) {
        (set<col>(row, values[col]), ...);
    }
    void set_row(size_t row, const std::array<ulint, num_cols>& values) {
        set_row(row, values, std::make_index_sequence<num_cols>{});
    }

    template<size_t... col>
    std::array<ulint, num_cols> get_row(size_t row, std::index_sequence<col...>) const { 
        return {get<col>(row)...};
    }
    std::array<ulint, num_cols> get_row(size_t row) const { 
        return get_row(row, std::make_index_sequence<num_cols>{});
    }


    [[nodiscard]] size_t size() const noexcept { return num_rows; }
    [[nodiscard]] size_t rows() const noexcept { return num_rows; }
    [[nodiscard]] static constexpr size_t cols() noexcept { return num_cols; }
    [[nodiscard]] size_t data_size() const noexcept {
        return byte_row_width * num_rows + sizeof(ulint) / sizeof(word_t);
    }
    [[nodiscard]] const std::array<uchar, num_cols>& get_widths() const noexcept { return widths; }

    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        out.write((char *)&num_rows, sizeof(num_rows));
        written_bytes += sizeof(num_rows);

        // Serialize column widths (may be zero columns)
        size_t widths_bytes = widths.size() * sizeof(uchar);
        if (widths_bytes > 0) {
            out.write((char *)widths.data(), widths_bytes);
            written_bytes += widths_bytes;
        }

        // Serialize packed data buffer
        const char* data_mem = reinterpret_cast<const char*>(data.data());
        size_t data_bytes = data.size() * sizeof(word_t);
        if (data_bytes > 0) {
            out.write(data_mem, data_bytes);
            written_bytes += data_bytes;
        }

        return written_bytes;
    }

    void load(std::istream &in) {
        in.read((char *)&num_rows, sizeof(num_rows));
        // Load column widths (may be zero columns)
        size_t widths_bytes = widths.size() * sizeof(uchar);
        if (widths_bytes > 0) {
            in.read((char *)widths.data(), widths_bytes);
        }
        init();
        // Load packed data buffer
        size_t data_bytes = data_size() * sizeof(word_t);
        if (data_bytes > 0) {
            in.read((char *)data.data(), data_bytes);
        }
    }

private:
    size_t num_rows;
    size_t byte_vector_width; // byte width of stored data (actual data size might be larger due to padding)
    size_t byte_row_width; // width of each row in bytes

    std::array<uchar, num_cols> widths; // Intended bit width of each column
    std::array<uchar, num_cols> byte_widths; // Actual byte width of each column
    std::array<uint16_t, num_cols> byte_offsets; // Offset of the first byte of each column in the vector
    std::array<ulint, num_cols> masks_extract; // Mask of the bits of each column for get
    std::array<ulint, num_cols> masks_write; // Mask of the bits of each column for set

    std::vector<word_t> data;

    inline size_t get_row_start_byte(size_t row) const {
        return row*byte_row_width;
    }

    void init() {
        size_t byte_pos = 0;
        for (size_t i = 0; i < num_cols; i++) {
            assert(widths[i] <= max_width);
            // Round up to the nearest byte
            if (widths[i] % num_bits_type(word_t) != 0) {
                byte_widths[i] = (widths[i] / num_bits_type(word_t)) + 1;
            } else {
                byte_widths[i] = widths[i] / num_bits_type(word_t);
            }

            masks_extract[i] = mask(widths[i]);
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

template<class columns>
class packed_vector_aligned : public packed_matrix_aligned<static_cast<size_t>(columns::COUNT)> {
    using base = packed_matrix_aligned<static_cast<size_t>(columns::COUNT)>;

public:
    packed_vector_aligned() = default;
    packed_vector_aligned(size_t rows, const std::array<uchar, static_cast<size_t>(columns::COUNT)>& widths)
        : base(rows, widths) {}

    template<columns col>
    ulint get(size_t row) const {
        return base::template get<static_cast<size_t>(col)>(row);
    }

    template<columns col>
    void set(size_t row, ulint val) {
        base::template set<static_cast<size_t>(col)>(row, val);
    }
};

class int_vector_aligned : public packed_matrix_aligned<1> {
    using base = packed_matrix_aligned<1>;
public:
    using value_type = ulint;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;

    class reference {
    public:
        reference(int_vector_aligned* v, size_type i) : vec(v), idx(i) {}

        reference(const reference&) = default;
        reference& operator=(const reference& other) {
            return *this = static_cast<ulint>(other);
        }

        reference& operator=(ulint value) {
            vec->set(idx, value);
            return *this;
        }

        operator ulint() const {
            return vec->get(idx);
        }

    private:
        int_vector_aligned* vec;
        size_type idx;
    };

    class const_reference {
    public:
        const_reference(const int_vector_aligned* v, size_type i) : vec(v), idx(i) {}

        operator ulint() const {
            return vec->get(idx);
        }

    private:
        const int_vector_aligned* vec;
        size_type idx;
    };

    class iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = int_vector_aligned::value_type;
        using difference_type   = int_vector_aligned::difference_type;
        using reference         = int_vector_aligned::reference;
        using pointer           = void;

        iterator() : vec(nullptr), idx(0) {}
        iterator(int_vector_aligned* v, size_type i) : vec(v), idx(i) {}

        reference operator*() const { return reference(vec, idx); }

        iterator& operator++() { ++idx; return *this; }
        iterator operator++(int) { iterator tmp(*this); ++(*this); return tmp; }

        iterator& operator--() { --idx; return *this; }
        iterator operator--(int) { iterator tmp(*this); --(*this); return tmp; }

        iterator& operator+=(difference_type n) { idx += n; return *this; }
        iterator& operator-=(difference_type n) { idx -= n; return *this; }

        iterator operator+(difference_type n) const { return iterator(vec, idx + n); }
        iterator operator-(difference_type n) const { return iterator(vec, idx - n); }

        difference_type operator-(const iterator& other) const {
            return static_cast<difference_type>(idx) - static_cast<difference_type>(other.idx);
        }

        bool operator==(const iterator& other) const { return vec == other.vec && idx == other.idx; }
        bool operator!=(const iterator& other) const { return !(*this == other); }
        bool operator<(const iterator& other) const { return idx < other.idx; }
        bool operator>(const iterator& other) const { return other < *this; }
        bool operator<=(const iterator& other) const { return !(other < *this); }
        bool operator>=(const iterator& other) const { return !(*this < other); }

        int_vector_aligned* vec;
        size_type idx;
    };

    class const_iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = int_vector_aligned::value_type;
        using difference_type   = int_vector_aligned::difference_type;
        using reference         = int_vector_aligned::const_reference;
        using pointer           = void;

        const_iterator() : vec(nullptr), idx(0) {}
        const_iterator(const int_vector_aligned* v, size_type i) : vec(v), idx(i) {}
        const_iterator(const iterator& it) : vec(it.vec), idx(it.idx) {}

        reference operator*() const { return reference(vec, idx); }

        const_iterator& operator++() { ++idx; return *this; }
        const_iterator operator++(int) { const_iterator tmp(*this); ++(*this); return tmp; }

        const_iterator& operator--() { --idx; return *this; }
        const_iterator operator--(int) { const_iterator tmp(*this); --(*this); return tmp; }

        const_iterator& operator+=(difference_type n) { idx += n; return *this; }
        const_iterator& operator-=(difference_type n) { idx -= n; return *this; }

        const_iterator operator+(difference_type n) const { return const_iterator(vec, idx + n); }
        const_iterator operator-(difference_type n) const { return const_iterator(vec, idx - n); }

        difference_type operator-(const const_iterator& other) const {
            return static_cast<difference_type>(idx) - static_cast<difference_type>(other.idx);
        }

        bool operator==(const const_iterator& other) const { return vec == other.vec && idx == other.idx; }
        bool operator!=(const const_iterator& other) const { return !(*this == other); }
        bool operator<(const const_iterator& other) const { return idx < other.idx; }
        bool operator>(const const_iterator& other) const { return other < *this; }
        bool operator<=(const const_iterator& other) const { return !(other < *this); }
        bool operator>=(const const_iterator& other) const { return !(*this < other); }

        const int_vector_aligned* vec;
        size_type idx;
    };

    int_vector_aligned() = default;
    int_vector_aligned(size_t rows, uchar width) : base(rows, {width}) {}
    int_vector_aligned(std::vector<ulint> data) 
    : base(data.size(), {(data.empty() ? static_cast<uchar>(0) : bit_width(*std::max_element(data.begin(), data.end())))}) {
        for (size_t i = 0; i < data.size(); i++) {
            set(i, data[i]);
        }
    }
    int_vector_aligned(std::vector<ulint> data, uchar width)
    : base(data.size(), {width}) {
        for (size_t i = 0; i < data.size(); i++) {
            set(i, data[i]);
        }
    }

    ulint get(size_t row) const {
        return base::template get<0>(row);
    }

    void set(size_t row, ulint val) {
        base::template set<0>(row, val);
    }

    reference operator[](size_t row) {
        return reference(this, row);
    }

    ulint operator[](size_t row) const {
        return get(row);
    }

    iterator begin() {
        return iterator(this, 0);
    }

    iterator end() {
        return iterator(this, this->size());
    }

    const_iterator begin() const {
        return cbegin();
    }

    const_iterator end() const {
        return cend();
    }

    const_iterator cbegin() const {
        return const_iterator(this, 0);
    }

    const_iterator cend() const {
        return const_iterator(this, this->size());
    }
};

inline void swap(int_vector_aligned::reference&& x, int_vector_aligned::reference&& y) noexcept {
    ulint t = static_cast<ulint>(x);
    x = static_cast<ulint>(y);
    y = t;
}

} // namespace orbit

#endif // end of include guard: _PACKED_VECTOR_ALIGNED_HPP