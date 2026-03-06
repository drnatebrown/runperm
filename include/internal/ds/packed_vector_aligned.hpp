#ifndef _PACKED_VECTOR_ALIGNED_HPP
#define _PACKED_VECTOR_ALIGNED_HPP

#include "internal/common.hpp"
#include <cassert>
#include <cstring>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <array>
#include <iostream>
#include <algorithm>
#include <iterator>

// Column values are rounded to the nearest byte
// TODO this would be better if definitely fast, using SIMD to unpack entire rows?
template <size_t NumCols>
class PackedMatrixAligned {
public:
    using word_t = uchar;

    // We read ulint at a time, this ensures we never need to read more than one ulint
    constexpr static uchar max_width = NUM_BITS(ulint); // should be 64 bits

    PackedMatrixAligned() = default;
    PackedMatrixAligned(const ulint rows, const std::array<uchar, NumCols>& widths) {
        PackedMatrixAligned::num_rows = rows;
        PackedMatrixAligned::widths = widths;

        init();
    }

    PackedMatrixAligned(PackedMatrixAligned&& other) noexcept = default;
    PackedMatrixAligned& operator=(PackedMatrixAligned&& other) noexcept = default;
    PackedMatrixAligned(const PackedMatrixAligned& other) = default;
    PackedMatrixAligned& operator=(const PackedMatrixAligned& other) = default;
    ~PackedMatrixAligned() = default;

    template<size_t Col>
    ulint get(size_t row) const {
        static_assert(Col < NumCols, "Column out of bounds");
        assert(row < num_rows);

        size_t byte_pos = get_row_start_byte(row) + byte_offsets[Col];
        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        return extract_bits(bits, masks_extract[Col]);
    }      

    template<size_t Col>
    void set(size_t row, ulint val) {
        static_assert(Col < NumCols, "Column out of bounds");
        assert(row < num_rows);
        assert(val < POW2(BYTES_TO_BITS(byte_widths[Col]))); // value must fit in the column

        size_t byte_pos = get_row_start_byte(row) + byte_offsets[Col];

        ulint bits = 0;
        std::memcpy(&bits, &data[byte_pos], sizeof(ulint));
        write_bits(bits, masks_write[Col], val);
        std::memcpy(&data[byte_pos], &bits, sizeof(ulint));
    }

    template<size_t... Col>
    void set_row(size_t row, const std::array<ulint, NumCols>& values, std::index_sequence<Col...>) {
        (set<Col>(row, values[Col]), ...);
    }
    void set_row(size_t row, const std::array<ulint, NumCols>& values) {
        set_row(row, values, std::make_index_sequence<NumCols>{});
    }

    template<size_t... Col>
    std::array<ulint, NumCols> get_row(size_t row, std::index_sequence<Col...>) const { 
        return {get<Col>(row)...};
    }
    std::array<ulint, NumCols> get_row(size_t row) const { 
        return get_row(row, std::make_index_sequence<NumCols>{});
    }


    [[nodiscard]] size_t size() const noexcept { return num_rows; }
    [[nodiscard]] size_t rows() const noexcept { return num_rows; }
    [[nodiscard]] static constexpr size_t cols() noexcept { return NumCols; }
    [[nodiscard]] size_t data_size() const noexcept {
        return byte_row_width * num_rows + sizeof(ulint) / sizeof(word_t);
    }
    [[nodiscard]] const std::array<uchar, NumCols>& get_widths() const noexcept { return widths; }

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

    std::array<uchar, NumCols> widths; // Intended bit width of each column
    std::array<uchar, NumCols> byte_widths; // Actual byte width of each column
    std::array<uint16_t, NumCols> byte_offsets; // Offset of the first byte of each column in the vector
    std::array<ulint, NumCols> masks_extract; // Mask of the bits of each column for get
    std::array<ulint, NumCols> masks_write; // Mask of the bits of each column for set

    std::vector<word_t> data;

    inline size_t get_row_start_byte(size_t row) const {
        return row*byte_row_width;
    }

    void init() {
        size_t byte_pos = 0;
        for (size_t i = 0; i < NumCols; i++) {
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

template<class Columns>
class PackedVectorAligned : public PackedMatrixAligned<static_cast<size_t>(Columns::COUNT)> {
    using Base = PackedMatrixAligned<static_cast<size_t>(Columns::COUNT)>;

public:
    PackedVectorAligned() = default;
    PackedVectorAligned(size_t rows, const std::array<uchar, static_cast<size_t>(Columns::COUNT)>& widths)
        : Base(rows, widths) {}

    template<Columns Col>
    ulint get(size_t row) const {
        return Base::template get<static_cast<size_t>(Col)>(row);
    }

    template<Columns Col>
    void set(size_t row, ulint val) {
        Base::template set<static_cast<size_t>(Col)>(row, val);
    }
};

class IntVectorAligned : public PackedMatrixAligned<1> {
    using Base = PackedMatrixAligned<1>;
public:
    using value_type = ulint;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;

    class reference {
    public:
        reference(IntVectorAligned* v, size_type i) : vec(v), idx(i) {}

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
        IntVectorAligned* vec;
        size_type idx;
    };

    class const_reference {
    public:
        const_reference(const IntVectorAligned* v, size_type i) : vec(v), idx(i) {}

        operator ulint() const {
            return vec->get(idx);
        }

    private:
        const IntVectorAligned* vec;
        size_type idx;
    };

    class iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = IntVectorAligned::value_type;
        using difference_type   = IntVectorAligned::difference_type;
        using reference         = IntVectorAligned::reference;
        using pointer           = void;

        iterator() : vec(nullptr), idx(0) {}
        iterator(IntVectorAligned* v, size_type i) : vec(v), idx(i) {}

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

        IntVectorAligned* vec;
        size_type idx;
    };

    class const_iterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = IntVectorAligned::value_type;
        using difference_type   = IntVectorAligned::difference_type;
        using reference         = IntVectorAligned::const_reference;
        using pointer           = void;

        const_iterator() : vec(nullptr), idx(0) {}
        const_iterator(const IntVectorAligned* v, size_type i) : vec(v), idx(i) {}
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

        const IntVectorAligned* vec;
        size_type idx;
    };

    IntVectorAligned() = default;
    IntVectorAligned(size_t rows, uchar width) : Base(rows, {width}) {}
    IntVectorAligned(std::vector<ulint> data) 
    : Base(data.size(), {(data.empty() ? static_cast<uchar>(0) : bit_width(*std::max_element(data.begin(), data.end())))}) {
        for (size_t i = 0; i < data.size(); i++) {
            set(i, data[i]);
        }
    }

    ulint get(size_t row) const {
        return Base::template get<0>(row);
    }

    void set(size_t row, ulint val) {
        Base::template set<0>(row, val);
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

// Enable std algorithms that rely on std::iter_swap for IntVectorAligned iterators.
namespace std {
    template<>
    inline void iter_swap(IntVectorAligned::iterator a, IntVectorAligned::iterator b) {
        ulint tmp = static_cast<ulint>(*a);
        *a = static_cast<ulint>(*b);
        *b = tmp;
    }
}

#endif // end of include guard: _PACKED_VECTOR_ALIGNED_HPP