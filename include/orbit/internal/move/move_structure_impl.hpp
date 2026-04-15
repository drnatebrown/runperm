#ifndef _INTERNAL_MOVE_STRUCTURE_HPP
#define _INTERNAL_MOVE_STRUCTURE_HPP

#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <numeric>

#include "orbit/common.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/move/move_table.hpp"
#include "orbit/internal/move/move_splitting.hpp"

namespace orbit {

template <typename columns_t = move_columns, template<typename> class table_t = move_vector>
class move_structure
{
public:
    // Sets NumCols, Columns, and cols_traits
    MOVE_CLASS_TRAITS(typename table_t<columns_t>::columns)
    using position = typename cols_traits::position;

    move_structure() = default;
    
    // Constructor from permutation data
    move_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& images, const split_params& split_params = split_params())
    : move_structure(interval_encoding_impl<>::from_lengths_and_images(lengths, images, split_params)) {}

    template<typename interval_encoding_t>
    move_structure(const interval_encoding_t& enc) {
        n = enc.domain();
        r = enc.runs();
        table = find_structure(enc);
    }

    // Constructor from pre-computed table (move semantics) for advanced users
    move_structure(packed_vector<columns> &&structure, const size_t domain, const ulint runs) 
        : table(std::move(structure)), n(domain), r(runs) {}

    template<typename interval_encoding_t>
    static packed_vector<columns> find_structure(const interval_encoding_t& enc) {
        packed_vector<columns> structure(enc.intervals(), get_move_widths(enc.domain(), enc.intervals(), enc.max_length()));
        populate_structure(structure, enc);
        return structure;
    }

    // === Interval Start/Length Accessors ===
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, ulint>
    get_start(size_t i) const {
        assert(i <= table.size());
        return (i == table.size()) ? n : table.get_start(i);
    }
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, ulint>
    inline get_start(position pos) const {
        return get_start(pos.interval);
    }

    ulint get_length(size_t i) const {
        assert(i < table.size());
        if constexpr (cols_traits::RELATIVE) {
            return table.get_length(i);
        } else {
            return (i == table.size() - 1)
                ? n - get_start(i)
                : get_start(i + 1) - get_start(i);
        }
    }
    inline ulint get_length(position pos) const {
        return get_length(pos.interval);
    }

    // === Pointer/Offset Accessors ===
    ulint get_pointer(size_t i) const {
        assert(i < table.size());
        return table.get_pointer(i);
    }
    inline ulint get_pointer(position pos) const {
        return get_pointer(pos.interval);
    }
    
    ulint get_offset(size_t i) const {
        assert(i < table.size());
        return table.get_offset(i);
    }
    inline ulint get_offset(position pos) const {
        return get_offset(pos.interval);
    }

    // === Generic Column Accessors ===
    template<columns col>
    ulint get(size_t i) const {
        return table.template get<col>(i);
    }
    template<columns col>
    ulint get(position pos) const {
        return get<col>(pos.interval);
    }
    std::array<ulint, num_cols> get_row(size_t i) const {
        return table.get_row(i);
    }
    std::array<ulint, num_cols> get_row(position pos) const {
        return get_row(pos.interval);
    }

    // === Structure Properties ===
    ulint domain() const {
        return n;
    }
    ulint runs() const {
        return r;
    }
    ulint intervals() const {
        return table.size();
    }

    const std::array<uchar, num_cols>& get_widths() const {
        return table.get_widths();
    }

    // === position Navigation ===
    position first() const {
        return position();
    }
    position last() const {
        size_t interval = table.size() - 1;
        size_t offset = get_length(interval) - 1;
        if constexpr (cols_traits::RELATIVE) {
            return position{interval, offset};
        } else {
            return position{interval, offset, n - 1};
        }
    }

    position move(position pos) const
    {
        assert(pos.interval < table.size());
        if constexpr (cols_traits::RELATIVE) {
            assert(pos.offset < get_length(pos));
            pos = {get_pointer(pos), pos.offset + get_offset(pos)};
        } else {
            assert(pos.idx < get_start(pos.interval + 1));
            ulint next_interval = get_pointer(pos);
            ulint next_offset = get_offset(pos) + pos.offset;
            pos = {next_interval, next_offset, get_start(next_interval) + next_offset};
        }
        return fast_forward(pos);
    }

    position move_exponential(position pos) const {
        if constexpr (cols_traits::RELATIVE) {
            return move(pos);
        } else {
            assert(pos.interval < table.size());
            assert(pos.idx < get_start(pos.interval + 1));
            ulint next_interval = get_pointer(pos);
            ulint next_offset = get_offset(pos) + pos.offset;
            pos = {next_interval, next_offset, get_start(next_interval) + next_offset};
            return fast_forward_exponential(pos);
        }
    }

    // === Utility Methods ===
    std::string get_file_extension() const
    {
        return MOVE_STRUCTURE_EXTENSION;
    }

    void move_stats() const
    {
        std::cout << "Number of contigious permutations: r = " << r << std::endl;
        std::cout << "Number of move intervals: r' = " << table.size() << std::endl;
        std::cout << "Domain of permutation: n = " << n << std::endl;
        std::cout << "Rate n/r = " << double(n) / r << std::endl;
    }

    // === Serialization ===
    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        written_bytes += table.serialize(out);

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t size;

        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        table.load(in);
    }

protected:
    table_t<columns_t> table;
    ulint n;
    ulint r;
    
    static std::array<uchar, num_cols> get_move_widths(const ulint domain, const ulint intervals, const ulint max_length) {
        std::array<uchar, num_cols> widths = {0};
        for (size_t i = 0; i < num_cols; ++i) {
            widths[i] = bytes_to_bits(DEFAULT_BYTES);
        }

        if constexpr (cols_traits::RELATIVE) {
            widths[static_cast<size_t>(cols_traits::PRIMARY)] = bit_width(max_length);
        } else {
            widths[static_cast<size_t>(cols_traits::PRIMARY)] = bit_width(domain - 1);
        }

        widths[static_cast<size_t>(cols_traits::POINTER)] = bit_width(intervals - 1);
        widths[static_cast<size_t>(cols_traits::OFFSET)] = bit_width(max_length);

        return widths;
    }
    
    template<typename interval_encoding_t>
    static void populate_structure(packed_vector<columns>& structure, const interval_encoding_t& enc) {
        size_t start_val = 0;
        size_t output_start_val = 0;
        size_t img_rank_inv_idx = 0;
        for (size_t i = 0; i < enc.intervals(); ++i) {
            size_t length = enc.get_length(i);
            if constexpr (cols_traits::RELATIVE) {
                structure.template set<to_cols(cols_traits::PRIMARY)>(i, length);
            }
            else {
                structure.template set<to_cols(cols_traits::PRIMARY)>(i, start_val);
            }

            while (img_rank_inv_idx < enc.intervals() && output_start_val < start_val + length) {
                ulint img_rank_inv_val = enc.get_img_rank_inv(img_rank_inv_idx);
                structure.template set<to_cols(cols_traits::POINTER)>(img_rank_inv_val, i);
                structure.template set<to_cols(cols_traits::OFFSET)>(img_rank_inv_val, output_start_val - start_val);
                output_start_val += enc.get_length(img_rank_inv_val);
                ++img_rank_inv_idx;
            }
            start_val += length;
        }
    }

    // === position Navigation ===
    inline position fast_forward(position pos) const {
        if constexpr (cols_traits::RELATIVE) {
            ulint length = get_length(pos);
            while (pos.offset >= length) {
                pos.offset -= length;
                ++pos.interval;
                length = get_length(pos);
            }    
        } else {
            ulint curr_start = pos.idx - pos.offset;
            ulint next_start = get_start(pos.interval + 1);
            while (pos.idx >= next_start) {
                pos.offset -= next_start - curr_start;
                ++pos.interval;
                curr_start = next_start;
                next_start = get_start(pos.interval + 1);
            }
        }
        return pos;
    }

    inline position fast_forward_exponential(position pos) const {
        if constexpr (cols_traits::RELATIVE) {
            return fast_forward(pos);
        } else {
            ulint curr_start = pos.idx - pos.offset;
            if (pos.idx < get_start(pos.interval + 1)) return pos;
        
            // Exponential search to find upper bound
            ulint hi = 1;
            while (pos.interval + hi < table.size() && pos.idx >= get_start(pos.interval + hi)) {
                hi *= 2;
            }
            if (pos.interval + hi > table.size()) hi = table.size() - pos.interval;
        
            // Binary search in [hi/2, hi] for largest k with pos.idx >= get_start(pos.interval + k)
            ulint lo = hi / 2;
            while (lo + 1 < hi) {
                ulint mid = lo + (hi - lo) / 2;
                if (pos.idx >= get_start(pos.interval + mid)) lo = mid;
                else hi = mid;
            }
            ulint k = lo;
        
            ulint new_start = get_start(pos.interval + k);
            pos.offset -= (new_start - curr_start);
            pos.interval += k;
            return pos;
        }
    }
};

} // namespace orbit

#endif /* end of include guard: _INTERNAL_MOVE_STRUCTURE_HPP */