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
#include "ds/packed_vector.hpp"

struct PermutationStats {
    size_t split_num_rows; // number of intervals after capping intervals to max allowed length
    size_t permutation_size; // n
    ulint max_observed_length; // max length of an interval, after splitting
    std::optional<ulint> max_allowed_length; // max length of an interval, used for splitting

    PermutationStats(size_t num_split_rows, size_t permutation_size, ulint max_observed_length, std::optional<ulint> max_allowed) :
        split_num_rows(num_split_rows), permutation_size(permutation_size), max_observed_length(max_observed_length), max_allowed_length(max_allowed) {}

    PermutationStats(const std::vector<ulint> &lengths, const std::optional<ulint> max_allowed = std::nullopt) {
        split_num_rows = lengths.size();
        permutation_size = 0;
        max_observed_length = 0;
        max_allowed_length = max_allowed;

        if (max_allowed_length) { 
            for (ulint length : lengths) {
                ulint num_splits = length / *max_allowed_length;
                permutation_size += length;
                if (num_splits > 0) {
                    split_num_rows += (num_splits - 1);
                    if (length % *max_allowed_length != 0) split_num_rows++;
                    max_observed_length = *max_allowed_length;
                }
                else {
                    max_observed_length = std::max(max_observed_length, length);
                }
            }
        } else {
            for (ulint length : lengths) {
                max_observed_length = std::max(max_observed_length, length);
                permutation_size += length;
            }
        }
    }
};

template<typename RunData,
         bool IntegratedMoveStructure,
         bool SplitIntervals,
         bool StoreAbsolutePositions,
         typename BaseColumns,
         template<typename> class TableType,
         template<typename> class StructureType,
         template<typename> class PackedType>
class RunPerm;

template <typename Table = MoveTable<>>
class MoveStructure
{
public:
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(typename Table::Columns)
    using Position = typename ColsTraits::Position;

    MoveStructure() = default;
    MoveStructure(
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& interval_permutation,
        std::optional<ulint> max_allowed_length = std::nullopt
    ) : MoveStructure(lengths, interval_permutation, 
                      PermutationStats(lengths, max_allowed_length)) {}
    
    // Advanced constructor from pre-built PackedVector
    MoveStructure(PackedVector<Columns>& move_permutation, size_t n, size_t r) 
        : table(move_permutation), n(n), r(r) {}

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

    template<Columns Col>
    ulint get(size_t i) const {
        return table.template get<Col>(i);
    }
    template<Columns Col>
    ulint get(Position pos) const {
        return get<Col>(pos.interval);
    }

    ulint size() const {
        return n;
    }

    ulint runs() const {
        return r;
    }

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

    template<typename RunData,
         bool IntegratedMoveStructure,
         bool SplitIntervals,
         bool StoreAbsolutePositions,
         typename BaseColumns,
         template<typename> class TableType,
         template<typename> class StructureType,
         template<typename> class PackedType>
    friend class RunPerm;

    // For friend classes who may have duplicated stats computation already, or after explicitly computing them
    MoveStructure(
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& interval_permutation,
        const PermutationStats& stats
    ) {
        r = lengths.size();
        n = stats.permutation_size;
        
        PackedVector<Columns> structure = find_structure(lengths, interval_permutation, stats);
        table = Table(structure);
    }

    // Used by friend classes who also set other data types using the callback
    // For this to work, need to pass the correct widths
    template<typename OnRowSetCallback>
    static MoveStructure construct_with_extended_columns(
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& interval_permutation,
        const PermutationStats& stats,
        const std::array<uchar, NumCols>& column_widths,
        OnRowSetCallback on_row_set
    ) {
        MoveStructure ms;
        ms.r = lengths.size();
        ms.n = stats.permutation_size;
        
        PackedVector<Columns> structure = find_structure(
            lengths, interval_permutation, stats, column_widths, on_row_set
        );
        
        ms.table = Table(structure);
        return ms;
    }
    
    template<typename Func>
    static void for_each_split_interval(
        const std::vector<ulint>& lengths,
        const ulint max_split_length,
        Func func  // func(orig_interval_idx, length)
    ) {
        for (size_t i = 0; i < lengths.size(); ++i) {
            if (lengths[i] > max_split_length) {
                ulint remaining = lengths[i];
                while (remaining > 0) {
                    ulint chunk = std::min(remaining, max_split_length);
                    func(i, chunk);
                    remaining -= chunk;
                }
            }
            else {
                func(i, lengths[i]);
            }
        }
    }
    
    static std::vector<ulint> compute_split_interval_permutations(
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& interval_permutation,
        const ulint split_num_rows,
        const ulint max_split_length
    ) {
        std::vector<ulint> split_interval_permutations(split_num_rows);
        size_t tbl_idx = 0;
        size_t curr_i = 0;
        size_t i_offset = 0;
        for_each_split_interval(lengths, max_split_length, [&](size_t i, size_t length) {
            if (i != curr_i) {
                i_offset = 0;
                curr_i = i;
            }
            split_interval_permutations[tbl_idx++] = interval_permutation[i] + i_offset;
            i_offset += length;
        });
        assert(tbl_idx == split_interval_permutations.size());
        return split_interval_permutations;
    }

    static std::vector<size_t> compute_sorted_indices(const std::vector<ulint>& interval_perm) {
        std::vector<size_t> indices(interval_perm.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return interval_perm[a] < interval_perm[b]; });
        return indices;
    }
    
    template<typename OnRowSetCallback>
    static void populate_structure(
        PackedVector<Columns>& structure,
        const std::vector<ulint>& lengths,
        const std::vector<ulint>& final_interval_permutation,
        const std::vector<size_t>& sorted_indices,
        const PermutationStats& stats,
        OnRowSetCallback on_row_set
    ) {
        size_t tbl_idx = 0;
        size_t start_val = 0;
        auto sort_itr = sorted_indices.begin();
        auto add_interval = [&](size_t i, size_t length) {
            if constexpr (ColsTraits::RELATIVE) {
                structure.template set<to_cols(ColsTraits::PRIMARY)>(tbl_idx, length);
            }
            else {
                structure.template set<to_cols(ColsTraits::PRIMARY)>(tbl_idx, start_val);
            }

            if constexpr (!std::is_same_v<OnRowSetCallback, std::nullptr_t>) {
                /* Callback with: structure (for setting extra columns), original interval row index (before splitting), 
                row index (in final structure, with splitting) */
                on_row_set(structure, i, tbl_idx);
            }

            while (sort_itr != sorted_indices.end() && final_interval_permutation[*sort_itr] < start_val + length) {
                structure.template set<to_cols(ColsTraits::POINTER)>(*sort_itr, tbl_idx);
                structure.template set<to_cols(ColsTraits::OFFSET)>(*sort_itr, final_interval_permutation[*sort_itr] - start_val);
                ++sort_itr;
            }
            ++tbl_idx;
            start_val += length;
        };

        if (stats.max_allowed_length) {
            for_each_split_interval(lengths, *stats.max_allowed_length, add_interval);
        }
        else {
            for (size_t i = 0; i < lengths.size(); ++i) {
                add_interval(i, lengths[i]);
            }
        }
    }

    // Finds the structure of permutations for move from the lengths and interval permutation of interval starts
    // callback to assist with setting other data in the table that isn't directly related to the move structure
    // this is necessary for, example, RunPerm, which hides this advanced logic from the user
    template<typename OnRowSetCallback>
    static PackedVector<Columns> find_structure(const std::vector<ulint> &lengths, const std::vector<ulint> &interval_permutation, const PermutationStats stats, std::array<uchar, NumCols> column_widths, OnRowSetCallback on_row_set) {
        PackedVector<Columns> structure(stats.split_num_rows, column_widths);

        std::vector<ulint> split_interval_permutation; // if computed, see below
        if (stats.max_allowed_length) {
            split_interval_permutation = compute_split_interval_permutations(lengths, interval_permutation, stats.split_num_rows, *stats.max_allowed_length);
        }
        // the actual final interval permutation, which is either the original or the split one
        const std::vector<ulint> *final_interval_permutation = stats.max_allowed_length ? &split_interval_permutation : &interval_permutation;

        // Sort the interval idx by their interval permutation, used to find pointers/offsets
        auto sorted_indices = compute_sorted_indices(*final_interval_permutation);
        
        populate_structure(structure, lengths, *final_interval_permutation, sorted_indices, stats, on_row_set);
        return structure;
    }

    static std::array<uchar, NumCols> get_move_widths(const PermutationStats& stats) {
        std::array<uchar, NumCols> widths = {0};
        for (size_t i = 0; i < NumCols; ++i) {
            widths[i] = BYTES_TO_BITS(DEFAULT_BYTES);
        }

        if constexpr (ColsTraits::RELATIVE) {
            widths[static_cast<size_t>(ColsTraits::PRIMARY)] = bit_width(stats.max_observed_length);
        } else {
            widths[static_cast<size_t>(ColsTraits::PRIMARY)] = bit_width(stats.permutation_size);
        }

        widths[static_cast<size_t>(ColsTraits::POINTER)] = bit_width(stats.split_num_rows);
        widths[static_cast<size_t>(ColsTraits::OFFSET)] = bit_width(stats.max_observed_length);

        return widths;
    }

    static PackedVector<Columns> find_structure(const std::vector<ulint> &lengths, const std::vector<ulint> &interval_permutation, const PermutationStats stats) {
        return find_structure(lengths, interval_permutation, stats, get_move_widths(stats), nullptr);
    }
};

using MoveStructureTbl = MoveStructure<MoveTable<>>;
using MoveStructureTblIdx = MoveStructure<MoveTableIdx>;
using MoveStructureVec = MoveStructure<MoveVector<>>;
using MoveStructureVecIdx = MoveStructure<MoveVectorIdx>;

#endif /* end of include guard: _Move_HPP */