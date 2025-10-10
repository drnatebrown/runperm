#ifndef _Move_HH
#define _Move_HH

#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <vector>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <numeric>

#include "common.hpp"
#include "move/move_table.hpp"
#include "move/move_splitting.hpp"
#include "ds/packed_vector.hpp"

template <typename Table = MoveTable<>>
class MoveStructure
{
public:
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(typename Table::Columns)
    using Position = typename ColsTraits::Position;

    MoveStructure() = default;
    
    // Constructor from permutation data
    MoveStructure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams& split_params = SplitParams()) 
        : table(find_structure(lengths, interval_permutation, domain, split_params)), n(domain), r(table.size()) {}

    // Constructor from pre-computed table (move semantics) for advanced users
    MoveStructure(PackedVector<Columns> &&structure, const ulint domain) 
        : table(std::move(structure)), n(domain), r(table.size()) {}

    static PackedVector<Columns> find_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        // Determined the final lengths and interval permutations, which change if we use splitting
        const std::vector<ulint>* final_lengths = &lengths;
        const std::vector<ulint>* final_interval_permutation = &interval_permutation;
        ulint max_length = 0;
        
        // Split the lengths and interval permutations if needed, using length/balancing splitting
        SplitResult split_result;
        apply_splitting(final_lengths, final_interval_permutation, max_length, split_params, split_result);
        
        // Initialize the structure with the final lengths and interval permutations
        PackedVector<Columns> structure(final_lengths->size(), get_move_widths(domain, final_lengths->size(), max_length));
        
        // Sort the interval idx by their interval permutation, used to find pointers/offsets
        auto sorted_indices = compute_sorted_indices(*final_interval_permutation);
        
        // Fill the structure with move permutation data
        populate_structure(structure, *final_lengths, *final_interval_permutation, sorted_indices);
        return structure;
    }

    // === Interval Start/Length Accessors ===
    template <typename C = Columns>
    std::enable_if_t<!ColsTraitsFor<C>::RELATIVE, ulint>
    get_start(size_t i) const {
        assert(i <= table.size());
        return (i == table.size()) ? n : table.get_start(i);
    }
    template <typename C = Columns>
    std::enable_if_t<!ColsTraitsFor<C>::RELATIVE, ulint>
    inline get_start(Position pos) const {
        return get_start(pos.interval);
    }

    ulint get_length(size_t i) const {
        assert(i < table.size());
        if constexpr (ColsTraits::RELATIVE) {
            return table.get_length(i);
        } else {
            return (i == table.size() - 1)
                ? n - get_start(i)
                : get_start(i + 1) - get_start(i);
        }
    }
    inline ulint get_length(Position pos) const {
        return get_length(pos.interval);
    }

    // === Pointer/Offset Accessors ===
    ulint get_pointer(size_t i) const {
        assert(i < table.size());
        return table.get_pointer(i);
    }
    inline ulint get_pointer(Position pos) const {
        return get_pointer(pos.interval);
    }
    
    ulint get_offset(size_t i) const {
        assert(i < table.size());
        return table.get_offset(i);
    }
    inline ulint get_offset(Position pos) const {
        return get_offset(pos.interval);
    }

    // === Generic Column Accessors ===
    template<Columns Col>
    ulint get(size_t i) const {
        return table.template get<Col>(i);
    }
    template<Columns Col>
    ulint get(Position pos) const {
        return get<Col>(pos.interval);
    }

    // === Structure Properties ===
    ulint size() const {
        return n;
    }
    ulint runs() const {
        return r;
    }

    // === Position Navigation ===
    Position first() const {
        return Position();
    }
    Position last() const {
        size_t interval = table.size() - 1;
        size_t offset = get_length(interval) - 1;
        if constexpr (ColsTraits::RELATIVE) {
            return Position(interval, offset);
        } else {
            return Position(interval, offset, n - 1);
        }
    }

    Position move(Position pos) const
    {
        assert(pos.interval < table.size());
        if constexpr (ColsTraits::RELATIVE) {
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

    // === Utility Methods ===
    std::string get_file_extension() const
    {
        return MOVE_STRUCTURE_EXTENSION;
    }

    void move_stats() const
    {
        std::cout << "Number of contigious permutations: r = " << r << std::endl;
        std::cout << "Length of permutation: n = " << n << std::endl;
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

private:
    Table table;
    ulint n;
    ulint r;

    // === Structure Building Helpers ===
    static void apply_splitting(const std::vector<ulint>*& final_lengths, const std::vector<ulint>*& final_interval_permutation, ulint& max_length, const SplitParams& split_params, SplitResult& split_result) {
        auto set_by_split_result = [&](SplitResult& split_result) {
            final_lengths = &split_result.lengths;
            final_interval_permutation = &split_result.interval_permutations;
            max_length = split_result.max_length;
        };
        if (split_params.max_allowed_length) {
            split_by_max_allowed_length(*final_lengths, *final_interval_permutation, *split_params.max_allowed_length, split_result);
            set_by_split_result(split_result);
        }
        if (split_params.balancing_factor) {
            split_by_balancing_factor(*final_lengths, *final_interval_permutation, *split_params.balancing_factor, split_result);
            set_by_split_result(split_result);
        }
        if (!split_params.max_allowed_length && !split_params.balancing_factor) {
            max_length = *std::max_element(final_lengths->begin(), final_lengths->end());
        }
    }
    
    static std::array<uchar, NumCols> get_move_widths(const ulint domain, const ulint runs, const ulint max_length) {
        std::array<uchar, NumCols> widths = {0};
        for (size_t i = 0; i < NumCols; ++i) {
            widths[i] = BYTES_TO_BITS(DEFAULT_BYTES);
        }

        if constexpr (ColsTraits::RELATIVE) {
            widths[static_cast<size_t>(ColsTraits::PRIMARY)] = bit_width(max_length);
        } else {
            widths[static_cast<size_t>(ColsTraits::PRIMARY)] = bit_width(domain);
        }

        widths[static_cast<size_t>(ColsTraits::POINTER)] = bit_width(runs);
        widths[static_cast<size_t>(ColsTraits::OFFSET)] = bit_width(max_length);

        return widths;
    }

    static std::vector<size_t> compute_sorted_indices(const std::vector<ulint>& interval_perm) {
        std::vector<size_t> indices(interval_perm.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return interval_perm[a] < interval_perm[b]; });
        return indices;
    }
    
    static void populate_structure(
        PackedVector<Columns>& structure,
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& interval_permutation,
        const std::vector<size_t>& sorted_indices
    ) {
        size_t tbl_idx = 0;
        size_t start_val = 0;
        auto sort_itr = sorted_indices.begin();
        for (size_t i = 0; i < lengths.size(); ++i) {
            size_t length = lengths[i];
            if constexpr (ColsTraits::RELATIVE) {
                structure.template set<to_cols(ColsTraits::PRIMARY)>(tbl_idx, length);
            }
            else {
                structure.template set<to_cols(ColsTraits::PRIMARY)>(tbl_idx, start_val);
            }

            while (sort_itr != sorted_indices.end() && interval_permutation[*sort_itr] < start_val + length) {
                structure.template set<to_cols(ColsTraits::POINTER)>(*sort_itr, tbl_idx);
                structure.template set<to_cols(ColsTraits::OFFSET)>(*sort_itr, interval_permutation[*sort_itr] - start_val);
                ++sort_itr;
            }
            ++tbl_idx;
            start_val += length;
        }
    }

    // === Position Navigation ===
    inline Position fast_forward(Position pos) const {
        if constexpr (ColsTraits::RELATIVE) {
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
};

using MoveStructureTbl = MoveStructure<MoveTable<>>;
using MoveStructureTblIdx = MoveStructure<MoveTableIdx>;
using MoveStructureVec = MoveStructure<MoveVector<>>;
using MoveStructureVecIdx = MoveStructure<MoveVectorIdx>;

#endif /* end of include guard: _Move_HPP */