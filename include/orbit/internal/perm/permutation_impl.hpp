#ifndef _INTERNAL_PERMUTATION_HPP
#define _INTERNAL_PERMUTATION_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/move_splitting.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/move/move_structure_impl.hpp"
#include "orbit/internal/perm/data_columns.hpp"

namespace orbit {

constexpr bool DEFAULT_INTEGRATED_MOVE_STRUCTURE = false;
constexpr bool DEFAULT_STORE_ABSOLUTE_POSITIONS = false;
constexpr bool DEFAULT_EXPONENTIAL_SEARCH = false; // Whether to use exponential search for next(), only used if store_absolute_positions is true

/* ============================================= Advanced Implemenation ============================================= */

// If we're integrating the run data alongside the move structure, we don't need to store it separately
template <typename data_columns_t, bool integrated_move_structure>
struct separated_data_holder;
template <typename data_columns_t>
struct separated_data_holder<data_columns_t, false> { [[no_unique_address]] packed_vector<data_columns_t> data_cols; };
template <typename data_columns_t>
struct separated_data_holder<data_columns_t, true> { /* empty */ };

// TODO InversePermutation, which builds both the forward and inverse move structures if needed
template<typename data_columns_t = empty_data_columns, // Fields to be stored alongside the move structure representing a runny permutation
         bool integrated_move_structure = DEFAULT_INTEGRATED_MOVE_STRUCTURE, // Whether to pack the run data alongside the move structure
         bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS, // Whether to store absolute positions instead of interval/offset
         bool exponential_search = store_absolute_positions && DEFAULT_EXPONENTIAL_SEARCH, // Whether to use exponential search for next() by default, only used if store_absolute_positions is true
         typename base_columns_t = move_columns,
         template<typename, template<typename> class> class move_structure_t = move_structure,
         template<typename> class table_t = move_vector>
         // TODO need PackedType option?
class permutation_impl : separated_data_holder<data_columns_t, integrated_move_structure> {
protected:
    // Helpful constants for number of base (move permutation information) columns and run (additional data) columns
    static constexpr size_t num_run_cols = static_cast<size_t>(data_columns_t::COUNT);
    static constexpr size_t num_base_cols = static_cast<size_t>(base_columns_t::COUNT);

    // Switch the base columns to use the correct relative/absolute indexing if needed
    using base_columns = switch_columns<base_columns_t, store_absolute_positions>;
    // Use data_tupleColumns if integrating user data alongside the move structure, otherwise just use the switched columns
    // data_columns_wrapper extends the base_columns traits to include the run data columns
    using columns_t = std::conditional_t<integrated_move_structure, 
        typename data_columns_wrapper<data_columns_t, base_columns>::e, 
        base_columns>;

    // Sets num_cols, columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)

    // Base structure type is the move structure without run data
    using move_structure_base = move_structure_t<base_columns, table_t>;
    using move_structure_perm = move_structure_t<columns, table_t>;

public:
    using data_columns = data_columns_t;
    using data_tuple = columns_tuple<data_columns>;
    using position = typename move_structure_perm::position;

    // TODO use int_vector and container templates here

    static_assert(has_count_enumerator<data_columns>::value, "data_columns_t must have a COUNT enumerator");
    static_assert(!(!store_absolute_positions && exponential_search), "Exponential search is only supported with absolute positions");

    // check if we're using move_table
    static constexpr bool is_move_table_type() {
        return std::is_same_v<table_t<void>, move_table<void>>;
    }

    // At compile time, check if trying to integrate user data to MoveTable, which uses bitpacked structs
    static_assert(!(integrated_move_structure && is_move_table_type()),
    "Cannot integrate user data with MoveTable. Use move_vector or set integrated_move_structure=false."
    "MoveTable for integrated data requires specialized implementation!");

    permutation_impl() = default;

    // Gated constructors for when data columns are empty
    // Constructor from permutation vector
    template<class dc = data_columns_t,
         std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    permutation_impl(std::vector<ulint>& perm, const split_params& sp = split_params())
    : permutation_impl([&]{
            auto [lengths, images] = get_permutation_intervals(perm);
            return interval_encoding_impl<>::from_lengths_and_images(lengths, images, sp);
    }()) {}

    template<typename interval_encoding_t,
             class dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    permutation_impl(const interval_encoding_t& enc) 
    : permutation_impl(enc, std::vector<data_tuple>(enc.intervals())) {}
    
    // Constructor from lengths and interval permutation
    template<typename container1_t, typename container2_t,
             class dc = data_columns_t,
             std::enable_if_t<std::is_same_v<dc, empty_data_columns>, int> = 0>
    permutation_impl(const container1_t& lengths, const container2_t& images, const split_params& sp = split_params())
    : permutation_impl(lengths, images, sp, std::vector<data_tuple>(lengths.size())) {}

    // Intended constructor for manual splitting of run data
    template<typename interval_encoding_t>
    permutation_impl(const interval_encoding_t& enc, const std::vector<data_tuple> &run_data) {
        split_params_ = enc.get_split_params();
        packed_vector<base_columns> base_structure = move_structure_base::find_structure(enc);
        if (run_data.size() == enc.intervals()) {
            populate_structure(std::move(base_structure), run_data, enc.domain(), enc.runs());
        }
        else if (run_data.size() == enc.runs()) {
            throw std::invalid_argument("Run data size is same as number of runs, not intervals after splitting; avoid splitting, manually split run data, or use permutation copy split.");
        } else {
            throw std::invalid_argument("Run data size must be the same as the number of intervals (user defined splits) or number of runs (no splitting).");
        }
    }

    /**
     * lengths -> length of each interval which permutes contiguously
     * images -> permutation position of the first position of each interval
     * domain -> domain of the permutation, i.e. a permutatation over 1..n has domain n
     * run_data -> run data for each interval, the size of this vector should be the same as the number of intervals
     */
     // When user doesn't pass splitting params without providing permutation object input, use NO_SPLITTING
    template<typename container1_t, typename container2_t>
     permutation_impl(const container1_t &lengths, const container2_t &images, const std::vector<data_tuple> &run_data)
     : permutation_impl(lengths, images, NO_SPLITTING, run_data) {}

    // If splitting, copy the run data by default
    template<typename container1_t, typename container2_t>
    permutation_impl(const container1_t &lengths, const container2_t &images, const split_params &sp, const std::vector<data_tuple> &run_data)
        : permutation_impl(lengths, images, sp, run_data,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    // Path for constructor above, see below constructor for more details
    // This exists as a fast path in case no splitting is set and we do not need to call extend_run_data
    template<typename container1_t, typename container2_t>
    permutation_impl(const container1_t &lengths, const container2_t &images, const split_params &sp, const std::vector<data_tuple> &run_data, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(lengths.size() == images.size());
        split_params_ = sp;

        // Find the base structure (move structure without run data)
        auto enc = interval_encoding_impl<>::from_lengths_and_images(lengths, images, sp);
        size_t domain = enc.domain();
        size_t runs = enc.runs();
        packed_vector<base_columns> base_structure = move_structure_base::find_structure(enc);

        if (split_params_ == NO_SPLITTING) {
            populate_structure(std::move(base_structure), run_data, domain, runs);
        } else {
            std::vector<data_tuple> final_run_data = extend_run_data(lengths, base_structure, domain, get_run_cols_data);
            populate_structure(std::move(base_structure), final_run_data, domain, runs);
        }
    }

    // Advanced constructor for users who want to specify how to set the run data for each column type by passing a function:
    /** Function Signature:
     * data_tuple(ulint, ulint, ulint, ulint)
     * - the original interval index
     * - the original interval length
     * - the absolute offset of the first element in the new interval
     * - the new interval length
     * and returns the values to be stored in the new interval (run data) for that new row
     */
    template<typename container1_t, typename container2_t>
    permutation_impl(const container1_t &lengths, const container2_t &images, const split_params &sp, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(lengths.size() == images.size());
        split_params_ = sp;

        // Find the base structure (move structure without run data)
        auto enc = interval_encoding_impl<>::from_lengths_and_images(lengths, images, sp);
        size_t domain = enc.domain();
        size_t runs = enc.runs();
        packed_vector<base_columns> base_structure = move_structure_base::find_structure(enc);

        std::vector<data_tuple> final_run_data = extend_run_data(lengths, base_structure, domain, get_run_cols_data);
        populate_structure(std::move(base_structure), final_run_data, domain, runs);
    }

    // from pre-computed table (move semantics) for advanced users with integrated move structure
    static permutation_impl from_structure(packed_vector<columns> &&structure, const size_t domain, const size_t runs) {
        static_assert(integrated_move_structure, "Cannot construct permutation with pre-computed permutation structure if not integrating user data with move structure");
        return permutation_impl(move_structure_perm(std::move(structure), domain, runs));
    }

    // from pre-computed table (move semantics) for advanced users without integrated move structure
    static permutation_impl from_structure(packed_vector<base_columns> &&structure, std::vector<data_tuple> &run_data, const size_t domain, const size_t runs) {
        static_assert(!integrated_move_structure, "Cannot construct permutation with pre-computed permutation structure if integrating user data with move structure");
        auto result = permutation_impl(move_structure_perm(std::move(structure), domain, runs));
        result.populate_run_data(std::move(structure), run_data);
        return result;
    }

    static permutation_impl from_move_structure(move_structure_perm &&ms) {
        static_assert(integrated_move_structure, "Cannot construct permutation with pre-computed move structure if not integrating user data with move structure");
        return permutation_impl(std::move(ms));
    }

    static permutation_impl from_move_structure(move_structure_perm &&ms, std::vector<data_tuple> &run_data) {
        assert(run_data.size() == ms.size());
        static_assert(!integrated_move_structure, "Cannot construct permutation with pre-computed move structure if integrating user data with move structure");
        auto result = permutation_impl(std::move(ms));
        auto run_cols_widths = get_data_cols_widths(run_data);
        result.fill_separated_data(run_data, run_cols_widths);
        return result;
    }

    /*** NAVIGATION METHODS ***/
    
    position next(position position) { 
        if constexpr (store_absolute_positions && exponential_search) {
            return move_structure.move_exponential(position);
        } else {
            return move_structure.move(position);
        }
    }

    position next(position position, ulint steps) {
        for (ulint i = 0; i < steps; ++i) {
            if constexpr (store_absolute_positions && exponential_search) {
                position = move_structure.move_exponential(position);
            } else {
                position = move_structure.move(position);
            }
        }
        return position;
    }

    // Non-exponential search version of next()
    position next_linear(position position) {
        return move_structure.move(position);
    }

    position next_linear(position position, ulint steps) {
        for (ulint i = 0; i < steps; ++i) {
            position = move_structure.move(position);
        }
        return position;
    }

    // Exponential search version of next()
    position next_exponential(position position) {
        return move_structure.move_exponential(position);
    }

    position next_exponential(position position, ulint steps) {
        for (ulint i = 0; i < steps; ++i) {
            position = move_structure.move_exponential(position);
        }
        return position;
    }

    // Set position to interval above in underlying move structure, or circularly wrap to the bottom if already at top
    position up(position position) {
        if (position.interval == 0)
        {
            return last();
        }
        --position.interval;
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (store_absolute_positions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    // Set position to interval below in underlying move structure, or circularly wrap to the top if already at bottom
    position down(position position) {
        if (position.interval == move_structure.intervals() - 1)
        {
            return first();
        }
        ++position.interval;
        position.offset = 0;
        if constexpr (store_absolute_positions) {
            position.idx = move_structure.get_start(position.interval);
        }
        return position;
    }

    // Returns row/offset of largest idx before or at position run with matching run data value
    template<data_columns col,
             class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    std::optional<position> pred(position position, ulint val) {
        if (get<col>(position) == val) {
            return position;
        }

        while (get<col>(position) != val)
        {
            if (position.interval == 0) return std::nullopt;
            --position.interval;
        }
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (store_absolute_positions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    // Returns row/offset of smallest idx after or at position run with matching run data value
    template<data_columns col,
             class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    std::optional<position> succ(position position, ulint val) {
        if (get<col>(position) == val) {
            return position;
        }

        while (get<col>(position) != val)
        {
            if (position.interval == move_structure.runs() - 1) return std::nullopt;
            ++position.interval;
        }
        position.offset = 0;
        if constexpr (store_absolute_positions) {
            position.idx = move_structure.get_start(position.interval);
        }
        return position;
    }

    /*** PROPERTY METHODS ***/

    position first() { return move_structure.first(); }
    position last() { return move_structure.last(); }

    ulint domain() const { return move_structure.domain(); }
    ulint runs() const { return move_structure.runs(); }
    ulint intervals() const { return move_structure.intervals(); }

    const std::array<uchar, num_cols>& get_widths() const {
        return move_structure.get_widths();
    }

    /*** DATA COLUMNS ACCESS METHODS ***/

    template<data_columns col,
             class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    ulint get(ulint interval) const {
        if constexpr (integrated_move_structure) {
            return move_structure.template get<cols_traits::template data_column<col>()>(interval);
        } else {
            return this->data_cols.template get<col>(interval);
        }
    }
    template<data_columns col,
             class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    ulint get(position position) const {
        return get<col>(position.interval);
    }

    template<class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    data_tuple get_row(ulint interval) const {
        if constexpr (integrated_move_structure) {
            auto full_row = move_structure.get_row(interval);
            data_tuple run_row;
            for (size_t i = 0; i < num_run_cols; ++i) {
                run_row[i] = full_row[num_base_cols + i];
            }
            return run_row;
        } else {
            return this->data_cols.get_row(interval);
        }   
    }
    template<class dc = data_columns_t,
             std::enable_if_t<!std::is_same_v<dc, empty_data_columns>, int> = 0>
    data_tuple get_row(position position) const {
        return get_row(position.interval);
    }

    ulint get_length(ulint interval) const {
        return move_structure.get_length(interval);
    }
    ulint get_length(position position) const {
        return get_length(position.interval);
    }
    
    split_params get_split_params() const {
        return split_params_;
    }

    /*** SERIALIZATION METHODS ***/
    
    size_t serialize(std::ostream& os) {
        size_t written_bytes = 0;

        written_bytes += serialize_version(os);
        written_bytes += split_params_.serialize(os);

        written_bytes += move_structure.serialize(os);
        if constexpr (!integrated_move_structure) {
            written_bytes += this->data_cols.serialize(os);
        }
        return written_bytes;
    }

    void load(std::istream& is) {
        auto [serialized_major, serialized_minor, serialized_patch] = load_version(is);
        if (serialized_major != VERSION_MAJOR || serialized_minor != VERSION_MINOR || serialized_patch != VERSION_PATCH) {
            // TODO handle version mismatches
        }
        split_params_.load(is);
        move_structure.load(is);
        if constexpr (!integrated_move_structure) {
            this->data_cols.load(is);
        }
    }

protected:
    move_structure_perm move_structure;
    split_params split_params_;

    template<base_columns col>
    ulint get_base_column(size_t row) const {
        return move_structure.template get<to_cols(col)>(row);
    }
    template<base_columns col>
    ulint get_base_column(position position) const {
        return get_base_column<col>(position.interval);
    }
    
    template<class Container>
    std::vector<data_tuple> extend_run_data(const Container& lengths, const packed_vector<base_columns>& structure, const size_t domain, std::function<data_tuple(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        std::vector<data_tuple> final_run_data(structure.size());
        auto get_structure_length = [&](size_t row) {
            if constexpr (store_absolute_positions) {
                return (row == structure.size() - 1) 
                    ? domain - structure.template get<cols_traits::PRIMARY>(row)
                    : structure.template get<cols_traits::PRIMARY>(row + 1) - structure.template get<cols_traits::PRIMARY>(row);
            } else {
                return structure.template get<cols_traits::PRIMARY>(row);
            }
        };

        size_t orig_start = 0;
        size_t new_idx = 0;
        size_t new_start = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            ulint orig_interval = i;
            ulint orig_interval_length = lengths[i];

            while (new_idx < structure.size() && new_start < orig_start + orig_interval_length) {
                size_t new_offset_from_orig_start = new_start - orig_start;
                size_t new_length = get_structure_length(new_idx);
                final_run_data[new_idx] = get_run_cols_data(orig_interval, orig_interval_length, new_offset_from_orig_start, new_length);
                ++new_idx;
                new_start += new_length;
            } 
            orig_start += orig_interval_length;
        }
        assert(new_idx == final_run_data.size());

        return final_run_data;
    }

    std::array<uchar, num_run_cols> get_data_cols_widths(const std::vector<data_tuple> &run_data) {
        data_tuple max_value = {};
        std::array<uchar, num_run_cols> run_cols_widths = {};
        for (size_t i = 0; i < num_run_cols; ++i) {
            for (size_t j = 0; j < run_data.size(); ++j) {
                max_value[i] = std::max(max_value[i], run_data[j][i]);
            }
        }
        for (size_t i = 0; i < num_run_cols; ++i) {
            run_cols_widths[i] = bit_width(max_value[i]);
        }
        return run_cols_widths;
    }

    static std::array<uchar, num_cols> get_widths(const std::array<uchar, num_base_cols>& base_widths, const std::array<uchar, num_run_cols>& run_cols_widths) {
        std::array<uchar, num_cols> widths;
        for (size_t i = 0; i < num_base_cols; ++i) {
            widths[i] = base_widths[i];
        }
        for (size_t i = 0; i < num_run_cols; ++i) {
            widths[num_base_cols + i] = run_cols_widths[i];
        }
        return widths;
    }

    std::array<ulint, num_cols> get_row(const std::array<ulint, num_base_cols>& base_row, const data_tuple& run_row) {
        std::array<ulint, num_cols> row;
        for (size_t i = 0; i < num_base_cols; ++i) {
            row[i] = base_row[i];
        }
        for (size_t i = 0; i < num_run_cols; ++i) {
            row[num_base_cols + i] = run_row[i];
        }
        return row;
    }

    void fill_separated_data(const std::vector<data_tuple>& run_data, const std::array<uchar, num_run_cols>& run_cols_widths) {
        this->data_cols = packed_vector<data_columns>(run_data.size(), run_cols_widths);
        for (size_t i = 0; i < run_data.size(); ++i) {
            this->data_cols.set_row(i, run_data[i]);
        }
    }

    // Sets move structure and run data from the base structure and run data
    void populate_structure(packed_vector<base_columns>&& base_structure, const std::vector<data_tuple>& run_data, const size_t domain, const size_t runs) {
        auto run_cols_widths = get_data_cols_widths(run_data);
        if constexpr (integrated_move_structure) {
            auto base_widths = base_structure.get_widths();
            auto widths = get_widths(base_widths, run_cols_widths);

            packed_vector<columns> final_structure(base_structure.size(), widths);
            for (size_t i = 0; i < final_structure.size(); ++i) {
                auto base_row = base_structure.get_row(i);
                auto run_row = run_data[i];
                auto row = get_row(base_row, run_row);
                final_structure.set_row(i, row);
            }
            move_structure = move_structure_perm(std::move(final_structure), domain, runs);
        } else {
            move_structure = move_structure_perm(std::move(base_structure), domain, runs);
            fill_separated_data(run_data, run_cols_widths);
        }
    }
};

// A helper alias around permutation_impl without any run data, essentially just a move_structure
template<bool store_absolute_positions = DEFAULT_STORE_ABSOLUTE_POSITIONS, bool exponential_search = DEFAULT_EXPONENTIAL_SEARCH, typename base_columns = move_columns, template<typename, template<typename> class> class move_structure_t = move_structure, template<typename> class table_t = move_vector>
using move_permutation_impl = permutation_impl<empty_data_columns, false, store_absolute_positions, exponential_search, base_columns, move_structure_t, table_t>;

} // namespace orbit

#endif /* end of include guard: _INTERNAL_PERMUTATION_HPP */
