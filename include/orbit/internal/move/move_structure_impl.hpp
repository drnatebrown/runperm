#ifndef _INTERNAL_MOVE_STRUCTURE_HPP
#define _INTERNAL_MOVE_STRUCTURE_HPP

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
#include <algorithm>

#include "orbit/common.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/move/move_table.hpp"
#include "orbit/internal/move/move_splitting.hpp"

namespace orbit {

template<typename columns_t, template<typename> class table_t, typename derived>
class move_structure_base {
public:
    // Sets NumCols, Columns, and cols_traits
    MOVE_CLASS_TRAITS(typename table_t<columns_t>::columns)
    using position = typename cols_traits::position;
    using interval_encoding_t = interval_encoding_impl<cols_traits::INVERTIBLE>;

    move_structure_base() = default;
    
    // Constructor from permutation data
    template<typename container1_t, typename container2_t>
    move_structure_base(const container1_t& lengths, const container2_t& images, const split_params& sp = split_params())
    : move_structure_base(interval_encoding_t::from_lengths_and_images(lengths, images, sp)) {}

    template<typename interval_encoding_impl_t>
    move_structure_base(const interval_encoding_impl_t& enc)
    : move_structure_base(find_structure(enc), enc.domain(), enc.runs()) {}

    // Constructor from pre-computed table (move semantics) for advanced users
    move_structure_base(packed_vector<columns> &&structure, const size_t domain, const ulint runs) 
        : table(std::move(structure)), n(domain), r(runs) {}

    template<typename interval_encoding_impl_t>
    static packed_vector<columns> find_structure(const interval_encoding_impl_t& enc) {
        static_assert(interval_encoding_impl_t::invertible_tag == cols_traits::INVERTIBLE, "Invertible type mismatch");
        packed_vector<columns> structure(enc.intervals(), derived::get_move_widths(enc.domain(), enc.intervals(), enc.max_length()));
        derived::populate_structure(structure, enc);
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

    // === Utility Methods ===
    std::string get_file_extension() const
    {
        return derived::MOVE_STRUCTURE_EXTENSION;
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

        written_bytes += write_magic(out, derived::MAGIC);
        written_bytes += serialize_version(out);

        out.write((char *)&n, sizeof(n));
        written_bytes += sizeof(n);

        out.write((char *)&r, sizeof(r));
        written_bytes += sizeof(r);

        written_bytes += table.serialize(out);

        return written_bytes;
    }

    void load(std::istream &in)
    {
        check_magic(in, derived::MAGIC);
        auto [serialized_major, serialized_minor, serialized_patch] = load_version(in);
        if (serialized_major != VERSION_MAJOR || serialized_minor != VERSION_MINOR || serialized_patch != VERSION_PATCH) {
            // TODO handle version mismatches
        }

        in.read((char *)&n, sizeof(n));
        in.read((char *)&r, sizeof(r));

        table.load(in);
    }

protected:
    table_t<columns_t> table;
    ulint n;
    ulint r;

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

template<typename columns_t = move_columns, template<typename> class table_t = move_vector>
class move_structure_impl : public move_structure_base<columns_t, table_t, move_structure_impl<columns_t, table_t>> {
public:
    using base = move_structure_base<columns_t, table_t, move_structure_impl<columns_t, table_t>>;
    using base::base;

    MOVE_CLASS_TRAITS(typename table_t<columns_t>::columns)
    using position = typename cols_traits::position;

    static_assert(!base::cols_traits::INVERTIBLE);

    // === Pointer/Offset Accessors ===
    ulint get_pointer(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_pointer(i);
    }
    inline ulint get_pointer(position pos) const {
        return get_pointer(pos.interval);
    }
    
    ulint get_offset(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_offset(i);
    }
    inline ulint get_offset(position pos) const {
        return get_offset(pos.interval);
    }

    // === Position Navigation ===
    template<bool linear=true>
    position move(position pos) const {
        assert(pos.interval < this->table.size());
        if constexpr (cols_traits::RELATIVE) {
            assert(pos.offset < this->get_length(pos));
            pos = {get_pointer(pos), pos.offset + get_offset(pos)};
        } else {
            assert(pos.idx < this->get_start(pos.interval + 1));
            ulint next_interval = get_pointer(pos);
            ulint next_offset = get_offset(pos) + pos.offset;
            pos = {next_interval, next_offset, this->get_start(next_interval) + next_offset};
        }
        if constexpr (linear) {
            return this->fast_forward(pos);
        } else {
            return this->fast_forward_exponential(pos);
        }
    }

    position move_exponential(position pos) const {
        return move<false>(pos);
    }

    // === Required from Base Class ===
    // OrBit Move Structure
    static constexpr std::array<char, MAGIC_BYTES> MAGIC = {'O', 'B', 'M', 'S'};
    inline constexpr static const char MOVE_STRUCTURE_EXTENSION[] = ".move";

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

    template<typename interval_encoding_impl_t>
    static void populate_structure(packed_vector<columns>& structure, const interval_encoding_impl_t& enc) {
        static_assert(interval_encoding_impl_t::invertible_tag == false, "Invertible type mismatch");
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
};

template<typename columns_t = invertible_columns, template<typename> class table_t = move_vector>
class invertible_structure_impl : public move_structure_base<columns_t, table_t, invertible_structure_impl<columns_t, table_t>> {
public:
    using base = move_structure_base<columns_t, table_t, invertible_structure_impl<columns_t, table_t>>;
    using base::base;

    MOVE_CLASS_TRAITS(typename table_t<columns_t>::columns)
    using position = typename cols_traits::position;

    static_assert(base::cols_traits::INVERTIBLE);

    // === Pointer/Offset Accessors ===
    ulint get_pointer_fwd(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_pointer_fwd(i);
    }
    inline ulint get_pointer_fwd(position pos) const {
        return get_pointer_fwd(pos.interval);
    }

    ulint get_pointer_inv(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_pointer_inv(i);
    }
    inline ulint get_pointer_inv(position pos) const {
        return get_pointer_inv(pos.interval);
    }   

    bool get_fwd_interval(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_fwd_interval(i);
    }
    inline bool get_fwd_interval(position pos) const {
        return get_fwd_interval(pos.interval);
    }

    bool get_inv_interval(size_t i) const {
        assert(i < this->table.size());
        return this->table.get_inv_interval(i);
    }
    inline bool get_inv_interval(position pos) const {
        return get_inv_interval(pos.interval);
    }

    // === Position Navigation ===
    template<bool linear=true>
    position move_fwd(position pos) const { return move<true, linear>(pos); }
    template<bool linear=true>
    position move_inv(position pos) const { return move<false, linear>(pos); }
    position move_fwd_exponential(position pos) const { return move<true, false>(pos); }
    position move_inv_exponential(position pos) const { return move<false, false>(pos); }

    template<bool is_fwd=true, bool linear=true>
    position move(position pos) const {
        auto orig_interval_func = [this](size_t i) -> bool {
            if constexpr (is_fwd)
            return get_fwd_interval(i);
            else
            return get_inv_interval(i);
        };
        assert(orig_interval_func(0));
        auto get_pointer_func = [this](position pos) -> ulint {
            if constexpr (is_fwd)
            return get_pointer_fwd(pos);
            else
            return get_pointer_inv(pos);
        };
        
        ulint new_offset = pos.offset;
        while (!orig_interval_func(pos.interval)) {
            --pos.interval;
            new_offset += this->get_length(pos.interval);
        }

        if constexpr (cols_traits::RELATIVE) {
            pos = {get_pointer_func(pos), new_offset};
        } else {
            pos = {get_pointer_func(pos), new_offset, this->get_start(pos.interval) + new_offset};
        }
        if constexpr (linear) {
            return this->fast_forward(pos);
        } else {
            return this->fast_forward_exponential(pos);
        }
    }

    template<bool is_fwd=true>
    position move_exponential(position pos) const {
        return move<is_fwd, false>(pos);
    }

    // === Required from Base Class ===
    // OrBit Invertible Structure
    static constexpr std::array<char, MAGIC_BYTES> MAGIC = {'O', 'B', 'I', 'S'};
    inline constexpr static const char MOVE_STRUCTURE_EXTENSION[] = ".imove";

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

        widths[static_cast<size_t>(cols_traits::POINTER_FWD)] = bit_width(intervals - 1);
        widths[static_cast<size_t>(cols_traits::POINTER_INV)] = bit_width(intervals - 1);

        widths[static_cast<size_t>(cols_traits::FWD_INTERVAL)] = 1;
        widths[static_cast<size_t>(cols_traits::INV_INTERVAL)] = 1;

        return widths;
    }

    template<typename interval_encoding_impl_t>
    static void populate_structure(packed_vector<columns>& structure, const interval_encoding_impl_t& enc) {
        static_assert(interval_encoding_impl_t::invertible_tag == true, "Invertible type mismatch");
        size_t input_start_val = 0;
        size_t output_start_val = 0;
        size_t img_rank_inv_idx = 0;
        for (size_t i = 0; i < enc.intervals(); ++i) {
            size_t length = enc.get_length(i);
            if constexpr (cols_traits::RELATIVE) {
                structure.template set<to_cols(cols_traits::PRIMARY)>(i, length);
            }
            else {
                structure.template set<to_cols(cols_traits::PRIMARY)>(i, input_start_val);
            }
            structure.template set<to_cols(cols_traits::FWD_INTERVAL)>(i, enc.get_is_fwd_interval(i));
            structure.template set<to_cols(cols_traits::INV_INTERVAL)>(i, enc.get_is_inv_interval(i));

            
            while (img_rank_inv_idx < enc.intervals() && output_start_val < input_start_val + length) {
                ulint img_rank_inv_val = enc.get_img_rank_inv(img_rank_inv_idx);
                if (enc.get_is_fwd_interval(img_rank_inv_val)) {
                    assert(output_start_val == input_start_val);
                    structure.template set<to_cols(cols_traits::POINTER_FWD)>(img_rank_inv_val, i);
                }
                if (enc.get_is_inv_interval(img_rank_inv_val)) {
                    assert(output_start_val == input_start_val);
                    structure.template set<to_cols(cols_traits::POINTER_INV)>(i, img_rank_inv_val);
                }
                output_start_val += enc.get_length(img_rank_inv_val);
                ++img_rank_inv_idx;
            }
            input_start_val += length;
        }
    }
};

template<typename columns_t = move_columns, template<typename> class table_t = move_vector>
using move_structure = std::conditional_t<resolve_cols_traits<columns_t>::type::INVERTIBLE, 
    invertible_structure_impl<columns_t, table_t>, 
    move_structure_impl<columns_t, table_t>>;

} // namespace orbit

#endif /* end of include guard: _INTERNAL_MOVE_STRUCTURE_HPP */