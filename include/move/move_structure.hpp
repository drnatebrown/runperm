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

struct SplitParams {
    std::optional<ulint> max_allowed_length;
    std::optional<ulint> balancing_factor;

    SplitParams() : max_allowed_length(std::nullopt), balancing_factor(std::nullopt) {}
};

struct LengthSplitResult {
    std::vector<ulint> lengths;
    std::vector<ulint> interval_permutations;
    ulint max_length;
};

LengthSplitResult split_by_max_allowed_length(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint max_allowed_length) {
    assert(lengths.size() == interval_permutation.size());

    LengthSplitResult result;
    result.lengths.reserve(static_cast<size_t>(1.5*lengths.size()));
    result.interval_permutations.reserve(static_cast<size_t>(1.5*interval_permutation.size()));

    std::vector<ulint> new_lengths;
    std::vector<ulint> new_interval_permutations;
    result.max_length = 0;

    for (size_t i = 0; i < lengths.size(); ++i) {
        if (lengths[i] > max_allowed_length) {
            ulint remaining = lengths[i];
            size_t sum_to_curr_chunk = 0;
            while (remaining > 0) {
                ulint chunk = std::min(remaining, max_allowed_length);
                new_lengths.push_back(chunk);
                new_interval_permutations.push_back(interval_permutation[i] + sum_to_curr_chunk);
                remaining -= chunk;
                sum_to_curr_chunk += chunk;
            }
            result.max_length = max_allowed_length;
        } else {
            new_lengths.push_back(lengths[i]);
            new_interval_permutations.push_back(interval_permutation[i]);
            result.max_length = std::max(result.max_length, lengths[i]);
        }
    }
    new_lengths.shrink_to_fit();
    new_interval_permutations.shrink_to_fit();
    return result;
}

template <typename Table = MoveTable<>>
class MoveStructure
{
public:
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(typename Table::Columns)
    using Position = typename ColsTraits::Position;

    MoveStructure() = default;
    MoveStructure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams& split_params = SplitParams()) 
        : table(find_structure(lengths, interval_permutation, domain, split_params)), n(domain), r(table.size()) {}

    // Advanced constructor for creating a move structure from a pre-computed table
    MoveStructure(PackedVector<Columns> &&structure, const ulint domain) 
        : table(std::move(structure)), n(domain), r(table.size()) {}

    static PackedVector<Columns> find_structure(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams& split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());

        ulint max_length = domain;
        
        const std::vector<ulint> *final_lengths = &lengths;
        const std::vector<ulint> *final_interval_permutation = &interval_permutation;

        std::vector<ulint> split_lengths;
        std::vector<ulint> split_interval_permutations;
        if (split_params.max_allowed_length) {
            auto length_split_result = split_by_max_allowed_length(lengths, interval_permutation, *split_params.max_allowed_length);
            split_lengths = std::move(length_split_result.lengths);
            split_interval_permutations = std::move(length_split_result.interval_permutations);
            max_length = length_split_result.max_length;

            final_lengths = &split_lengths;
            final_interval_permutation = &split_interval_permutations;
        }
        if (split_params.balancing_factor) {
            // TODO
        }

        if (!split_params.max_allowed_length || split_params.balancing_factor) {
            max_length = *std::max_element(final_lengths->begin(), final_lengths->end());
        }

        return find_structure(*final_lengths, *final_interval_permutation, get_move_widths(domain, final_lengths->size(), max_length));
    }

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

    // Finds the structure of permutations for move from the lengths and interval permutation of interval starts
    static PackedVector<Columns> find_structure(const std::vector<ulint> &lengths, const std::vector<ulint> &interval_permutation, std::array<uchar, NumCols> column_widths) {
        PackedVector<Columns> structure(lengths.size(), column_widths);

        // Sort the interval idx by their interval permutation, used to find pointers/offsets
        auto sorted_indices = compute_sorted_indices(interval_permutation);
        
        populate_structure(structure, lengths, interval_permutation, sorted_indices);
        return structure;
    }
};

using MoveStructureTbl = MoveStructure<MoveTable<>>;
using MoveStructureTblIdx = MoveStructure<MoveTableIdx>;
using MoveStructureVec = MoveStructure<MoveVector<>>;
using MoveStructureVecIdx = MoveStructure<MoveVectorIdx>;

#endif /* end of include guard: _Move_HPP */