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

#include "common.hpp"
#include "ds/move_table.hpp"

struct PermutationStats {
    size_t split_num_rows; // number of intervals after capping intervals to max allowed length
    size_t permutation_size; // n
    ulint max_observed_length; // max length of an interval, after splitting
    std::optional<ulint> max_allowed_length; // max length of an interval, used for splitting

    PermutationStats(size_t num_rows, size_t permutation_size, ulint max_observed_length, std::optional<ulint> max_allowed) :
        split_num_rows(num_rows), permutation_size(permutation_size), max_observed_length(max_observed_length), max_allowed_length(max_allowed) {}
    
    PermutationStats(const std::vector<ulint> &lengths, const std::optional<ulint> max_allowed) {
        split_num_rows = lengths.size();
        permutation_size = 0;
        max_observed_length = 0;
        max_allowed_length = max_allowed;

        if (max_allowed_length) { 
            for (ulint length : lengths) {
                ulint num_splits = length / *max_allowed_length;
                split_num_rows += num_splits;
                max_observed_length = (num_splits > 0) ? *max_allowed_length : std::max(max_observed_length, length);
                permutation_size += length;

            }
        } else {
            for (ulint length : lengths) {
                max_observed_length = std::max(max_observed_length, length);
                permutation_size += length;
            }
        }
    }
};

template <typename Table = MoveTable<>>
class MoveStructure
{
public:
    using Columns = typename Table::Columns;
    using Position = typename MoveColsTraits<Columns>::Position;

    MoveStructure() {}

    // n and r are the size of the permutation and the number of positions such that either i == 0 or π(i-1) != π(i) - 1
    // Should be with respect to the original permutation, not the split permutation
    MoveStructure(PackedVector<Columns> &move_permutation, size_t n, size_t r) 
    : n(n), r(r), table(move_permutation) {}

    // As below, but assumes permutation stats were precomputed
    MoveStructure(const std::vector<ulint> &lengths, const std::vector<ulint> &interval_permutation, PermutationStats stats) 
    : MoveStructure(find_structure(lengths, interval_permutation, stats), stats.permutation_size, lengths.size()) {}

    // Length of interval and the permutation of the first position in each interval. If max allowed length is set, split intervals if greater than max allowed length
    MoveStructure(const std::vector<ulint> &lengths, const std::vector<ulint> &interval_permutation, std::optional<ulint> max_allowed_length = std::nullopt) 
    : MoveStructure(lengths, interval_permutation, PermutationStats(lengths, max_allowed_length)) {}

    ulint get_length(size_t i) const {
        assert(i < table.size());
        if constexpr (MoveColsTraits<Columns>::IS_LENGTH) {
            return table.get_length(i);
        } else if constexpr (MoveColsTraits<Columns>::IS_START) {
            return (i == table.size() - 1) ? n - table.get_start(i) : table.get_start(i + 1) - table.get_start(i);
        }
    }
    inline ulint get_length(Position pos) const {
        return get_length(pos.interval);
    }

    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_START, ulint>
    get_start(size_t i) const {
        assert(i <= table.size());
        return (i == table.size()) ? n : table.get_start(i);
    }
    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_START, ulint>
    inline get_start(Position pos) const {
        return get_start(pos.interval);
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

    ulint size() const {
        return n;
    }

    ulint runs() const {
        return r;
    }

    Position move(Position pos) const
    {
        assert(pos.interval < table.size());
        if constexpr (MoveColsTraits<Columns>::IS_LENGTH) {
            assert(pos.offset < get_length(pos));
            pos = {get_pointer(pos), pos.offset + get_offset(pos)};
        } else if constexpr (MoveColsTraits<Columns>::IS_START) {
            assert(pos.idx < get_start(pos.interval + 1));
            ulint next_interval = get_pointer(pos);
            pos = {next_interval, get_start(next_interval) + (pos.idx - get_start(pos)) + get_offset(pos)};
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

    size_t serialize(std::ostream &out) const
    {
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

    // Finds the structure of permutations for move from the lengths and interval permutation of interval starts
    static PackedVector<Columns> find_structure(const std::vector<ulint> lengths, const std::vector<ulint> interval_permutation, const PermutationStats stats) {
        std::array<uchar, Columns::NUM_COLS> widths = {0};
        if constexpr (MoveColsTraits<Columns>::IS_LENGTH) {
            widths[static_cast<size_t>(Columns::LENGTH)] = bit_width(stats.max_observed_length);
        } else if constexpr (MoveColsTraits<Columns>::IS_START) {
            widths[static_cast<size_t>(Columns::START)] = bit_width(stats.permutation_size);
        }
        widths[static_cast<size_t>(Columns::POINTER)] = bit_width(stats.split_num_rows);
        widths[static_cast<size_t>(Columns::OFFSET)] = bit_width(stats.max_observed_length);
        PackedVector<Columns> structure(stats.split_num_rows, widths);

        std::vector<ulint> *final_interval_permutation = &interval_permutation;

        auto split_loop = [&](auto func) {
            for (size_t i = 0; i < lengths.size(); ++i) {
                ulint max_splits = lengths[i]/ *stats.max_allowed_length;
                for (size_t split = 0; split < max_splits; ++split) {
                    func(i, *stats.max_allowed_length, split);
                }
                if (lengths[i] % *stats.max_allowed_length != 0) {
                    func(i, lengths[i] % *stats.max_allowed_length, max_splits);
                }
            }
        };

        // If intervals are split, need to find where the split intervals permute to
        // TODO can we handle this when we preprocessed the lengths?
        if (stats.max_allowed_length) {
            std::vector<ulint> new_interval_permutations = std::vector<ulint>(stats.split_num_rows);
            size_t tbl_idx = 0;
            auto add_interval_permutation = [&](size_t i, size_t _, size_t split) {
                new_interval_permutations[tbl_idx++] = interval_permutation[i] + split * (*stats.max_allowed_length);
            };
            split_loop(add_interval_permutation);
            final_interval_permutation = &new_interval_permutations;
        }

        // Sort the interval idx by their interval permutation
        std::vector<size_t> indices(final_interval_permutation->size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) { return (*final_interval_permutation)[a] < (*final_interval_permutation)[b]; });
        
        size_t tbl_idx = 0;
        size_t sort_idx = 0;
        size_t start_val = 0;
        auto add_interval = [&](size_t length) {
            if constexpr (MoveColsTraits<Columns>::IS_LENGTH) {
                structure.set<Columns::LENGTH>(tbl_idx, length);
            }
            else if constexpr (MoveColsTraits<Columns>::IS_START) {
                structure.set<Columns::START>(tbl_idx, start_val);
            }

            while (sort_idx < indices.size() && (*final_interval_permutation)[indices[sort_idx]] < start_val + length) {
                structure.set<Columns::POINTER>(indices[sort_idx], tbl_idx);
                structure.set<Columns::OFFSET>(indices[sort_idx], (*final_interval_permutation)[indices[sort_idx]] - start_val);
                ++sort_idx;
            }
            
            ++tbl_idx;
            start_val += length;
        };
        auto add_interval_split = [&](size_t _, size_t length, size_t _) {
            add_interval(length);
        };
        if (stats.max_allowed_length) {
            split_loop(add_interval_split);
        }
        else {
            for (size_t i = 0; i < lengths.size(); ++i) {
                add_interval(lengths[i]);
            }
        }

        return structure;
    }
    
    inline Position fast_forward(Position pos) const {
        if constexpr (MoveColsTraits<Columns>::IS_LENGTH) {
            while (pos.offset >= get_length(pos)) {
                pos.offset -= get_length(pos.interval++);
            }    
        } else if constexpr (MoveColsTraits<Columns>::IS_START) {
            while (pos.idx <= get_start(pos.interval + 1)) {
                ++pos.interval;
            }
        }
        return pos;
    }
};

#endif /* end of include guard: _Move_HPP */