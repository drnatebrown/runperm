#ifndef _INTERNAL_RUNPERM_HPP
#define _INTERNAL_RUNPERM_HPP

#include "internal/common.hpp"
#include "internal/move/move_structure.hpp"
#include "internal/runperm/run_columns.hpp"

constexpr bool DEFAULT_INTEGRATED_MOVE_STRUCTURE = false;
constexpr bool DEFAULT_STORE_ABSOLUTE_POSITIONS = false;

/* ============================================= Advanced Implemenation ============================================= */

// If we're integrating the run data alongside the move structure, we don't need to store it separately
template <typename RunColsType, bool IntegratedMoveStructure>
struct SeperatedDataHolder;
template <typename RunColsType>
struct SeperatedDataHolder<RunColsType, false> { [[no_unique_address]] PackedVector<RunColsType> run_cols_data; };
template <typename RunColsType>
struct SeperatedDataHolder<RunColsType, true> { /* empty */ };

// TODO InversePermutation, which builds both the forward and inverse move structures if needed
template<typename RunColsType, // Fields to be stored alongside the move structure representing a runny permutation
         bool IntegratedMoveStructure = DEFAULT_INTEGRATED_MOVE_STRUCTURE, // Whether to pack the run data alongside the move structure
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, // Whether to store absolute positions instead of interval/offset
         typename BaseColumnsType = MoveCols,
         template<typename, template<typename> class> class MoveStructureType = MoveStructure,
         template<typename> class TableType = MoveVector>
         // TODO need PackedType option?
class RunPermImpl : SeperatedDataHolder<RunColsType, IntegratedMoveStructure> {
protected:
    // Helpful constants for number of base (move permutation information) columns and run (additional data) columns
    static constexpr size_t NumRunCols = static_cast<size_t>(RunColsType::COUNT);
    static constexpr size_t NumBaseCols = static_cast<size_t>(BaseColumnsType::COUNT);

    // Switch the base columns to use the correct relative/absolute indexing if needed
    using BaseColumns = SwitchColumns<BaseColumnsType, StoreAbsolutePositions>;
    // Use RunDataColumns if integrating user data alongside the move structure, otherwise just use the switched columns
    // RunColsWrapper extends the BaseColumns traits to include the run data columns
    using ColumnsType = std::conditional_t<IntegratedMoveStructure, 
        typename RunColsWrapper<RunColsType, BaseColumns>::E, 
        BaseColumns>;

    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(ColumnsType)

    // Base structure type is the move structure without run data
    using MoveStructureBase = MoveStructureType<BaseColumns, TableType>;
    using MoveStructurePerm = MoveStructureType<Columns, TableType>;

public:
    using RunCols = RunColsType;
    using RunData = std::array<ulint, NumRunCols>;
    using Position = typename MoveStructurePerm::Position;

    // check if we're using MoveTable
    static constexpr bool is_move_table_type() {
        return std::is_same_v<TableType<void>, MoveTable<void>>;
    }

    // At compile time, check if trying to integrate user data to MoveTable, which uses bitpacked structs
    static_assert(!(IntegratedMoveStructure && is_move_table_type()),
    "Cannot integrate user data with MoveTable. Use MoveVector or set IntegratedMoveStructure=false."
    "MoveTable for integrated data requires specialized implementation!");

    RunPermImpl() = default;

    /**
     * lengths -> length of each interval which permutes contiguously
     * interval_permutation -> permutation position of the first position of each interval
     * domain -> domain of the permutation, i.e. a permutatation over 1..n has domain n
     * run_data -> run data for each interval, the size of this vector should be the same as the number of intervals
     */
    RunPermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const std::vector<RunData> &run_data) {
        assert(lengths.size() == interval_permutation.size());
        assert(run_data.size() == lengths.size());
        orig_intervals = lengths.size();
        position = Position(); // should start at 0

        // Find the base structure (move structure without run data)
        PackedVector<BaseColumns> base_structure = MoveStructureBase::find_structure(lengths, interval_permutation, domain);
        populate_structure(std::move(base_structure), run_data, domain);
    }

    // When Splitting, by default just copy the run data for the original interval if the move structure intervals have been split
    RunPermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams &split_params, const std::vector<RunData> &run_data)
        : RunPermImpl(lengths, interval_permutation, domain, split_params,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    // Advanced constructor for users who want to specify how to set the run data for each column type by passing a function:
    /** Function Signature:
     * RunData(ulint, ulint, ulint, ulint)
     * - the original interval index
     * - the original interval length
     * - the absolute offset of the first element in the new interval
     * - the new interval length
     * and returns the values to be stored in the new interval (run data) for that new row
     */
    RunPermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const SplitParams &split_params, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(lengths.size() == interval_permutation.size());
        orig_intervals = lengths.size();
        position = Position(); // should start at 0
        
        // Find the base structure (move structure without run data)
        PackedVector<BaseColumns> base_structure = MoveStructureBase::find_structure(lengths, interval_permutation, domain, split_params);
        std::vector<RunData> final_run_data = extend_run_data(lengths, domain, base_structure, get_run_cols_data);
        populate_structure(std::move(base_structure), final_run_data, domain);
    }

    // Constructor from pre-computed table (move semantics) for advanced users with integrated move structure
    RunPermImpl(PackedVector<Columns> &&structure, const ulint domain) : move_structure(MoveStructurePerm(std::move(structure), domain)), position(Position()), orig_intervals(structure.size()) {
        static_assert(IntegratedMoveStructure, "Cannot construct RunPerm with pre-computed permutation structure if not integrating user data with move structure");
    }

    // Constructor from pre-computed table (move semantics) for advanced users without integrated move structure
    RunPermImpl(PackedVector<BaseColumns> &&structure, std::vector<RunData> &run_data, const ulint domain) : move_structure(MoveStructurePerm(std::move(structure), domain)), position(Position()), orig_intervals(structure.size()) {
        static_assert(!IntegratedMoveStructure, "Cannot construct RunPerm with pre-computed permutation structure if integrating user data with move structure");
        populate_run_data(std::move(structure), run_data, domain);
    }

    RunPermImpl(MoveStructurePerm &&ms, const ulint domain) : move_structure(std::move(ms)), position(Position()), orig_intervals(move_structure.size()) {
        static_assert(IntegratedMoveStructure, "Cannot construct RunPerm with pre-computed move structure if not integrating user data with move structure");
    }

    RunPermImpl(MoveStructurePerm &&ms, std::vector<RunData> &run_data, const ulint domain) : move_structure(std::move(ms)), position(Position()), orig_intervals(move_structure.size()) {
        static_assert(!IntegratedMoveStructure, "Cannot construct RunPerm with pre-computed move structure if integrating user data with move structure");
        auto run_cols_widths = get_run_cols_widths(run_data);
        fill_seperated_data(run_data, run_cols_widths);
    }
    
    void next() { position = move_structure.move(position); }
    void next(ulint steps) {
        for (ulint i = 0; i < steps; ++i) {
            next();
        }
    }

    // Set position to interval above in underlying move structure, returns false if already at top
    bool up() {
        if (position.interval == 0)
        {
            return false;
        }
        --position.interval;
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return true;
    }

    // Set position to interval below in underlying move structure, returns false if already at bottom
    bool down() {
        if (position.interval == move_structure.runs() - 1)
        {
            return false;
        }
        ++position.interval;
        position.offset = 0;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval);
        }
        return true;
    }

    // Returns row/offset of largest idx before or at position run with matching run data value
    template<RunCols Col>
    std::optional<Position> pred(ulint val) {
        while (get<Col>() != val) 
        {
            if (position.interval == 0) return std::nullopt;
            --position.interval;
        }
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    // Returns row/offset of smallest idx after or at position run with matching run data value
    template<RunCols Col>
    std::optional<Position> succ(ulint val) {
        while (get<Col>() != val) 
        {
            if (position.interval == move_structure.runs() - 1) return std::nullopt;
            ++position.interval;
        }
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    void first() { position = move_structure.first(); }
    void last() { position = move_structure.last(); }

    Position get_position() const { return position; }
    void set_position(Position pos) { position = pos; }

    ulint size() const { return move_structure.size(); }
    ulint move_runs() const { return move_structure.runs(); }
    ulint permutation_runs() const { return orig_intervals; }

    template<RunCols Col>
    ulint get() const {
        if constexpr (IntegratedMoveStructure) {
            return move_structure.template get<ColsTraits::template run_column<Col>()>(position.interval);
        } else {
            return this->run_cols_data.template get<Col>(position.interval);
        }
    }

    ulint get_length() const {
        return move_structure.get_length(position.interval);
    }

    size_t serialize(std::ostream& os) {
        size_t written_bytes = 0;

        os.write((char *)&orig_intervals, sizeof(orig_intervals));
        written_bytes += sizeof(orig_intervals);

        written_bytes += move_structure.serialize(os);
        if constexpr (!IntegratedMoveStructure) {
            written_bytes += this->run_cols_data.serialize(os);
        }
        return written_bytes;
    }

    void load(std::istream& is) {
        is.read((char *)&orig_intervals, sizeof(orig_intervals));

        move_structure.load(is);
        if constexpr (!IntegratedMoveStructure) {
            this->run_cols_data.load(is);
        }
    }

protected:
    MoveStructurePerm move_structure;
    Position position;
    size_t orig_intervals; // before splitting, underlying move structure may be larger

    template<BaseColumns Col>
    ulint get_base_column() const {
        return move_structure.template get<Col>(position.interval);
    }
    
    std::vector<RunData> extend_run_data(const std::vector<ulint>& lengths, const ulint domain, const PackedVector<BaseColumns>& structure, std::function<RunData(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        std::vector<RunData> final_run_data(structure.size());
        auto get_structure_length = [&](size_t row) {
            if constexpr (StoreAbsolutePositions) {
                return (row == structure.size() - 1) 
                    ? domain - structure.template get<ColsTraits::PRIMARY>(row)
                    : structure.template get<ColsTraits::PRIMARY>(row + 1) - structure.template get<ColsTraits::PRIMARY>(row);
            } else {
                return structure.template get<ColsTraits::PRIMARY>(row);
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

    std::array<uchar, NumRunCols> get_run_cols_widths(const std::vector<RunData> &run_data) {
        RunData max_value = {};
        std::array<uchar, NumRunCols> run_cols_widths = {};
        for (size_t i = 0; i < NumRunCols; ++i) {
            for (size_t j = 0; j < run_data.size(); ++j) {
                max_value[i] = std::max(max_value[i], run_data[j][i]);
            }
        }
        for (size_t i = 0; i < NumRunCols; ++i) {
            run_cols_widths[i] = bit_width(max_value[i]);
        }
        return run_cols_widths;
    }

    static std::array<uchar, NumCols> get_widths(const std::array<uchar, NumBaseCols>& base_widths, const std::array<uchar, NumRunCols>& run_cols_widths) {
        std::array<uchar, NumCols> widths;
        for (size_t i = 0; i < NumBaseCols; ++i) {
            widths[i] = base_widths[i];
        }
        for (size_t i = 0; i < NumRunCols; ++i) {
            widths[NumBaseCols + i] = run_cols_widths[i];
        }
        return widths;
    }

    std::array<ulint, NumCols> get_row(const std::array<ulint, NumBaseCols>& base_row, const RunData& run_row) {
        std::array<ulint, NumCols> row;
        for (size_t i = 0; i < NumBaseCols; ++i) {
            row[i] = base_row[i];
        }
        for (size_t i = 0; i < NumRunCols; ++i) {
            row[NumBaseCols + i] = run_row[i];
        }
        return row;
    }

    void fill_seperated_data(const std::vector<RunData>& run_data, const std::array<uchar, NumRunCols>& run_cols_widths) {
        this->run_cols_data = PackedVector<RunCols>(run_data.size(), run_cols_widths);
        for (size_t i = 0; i < run_data.size(); ++i) {
            this->run_cols_data.template set_row(i, run_data[i]);
        }
    }

    // Sets move structure and run data from the base structure and run data
    void populate_structure(PackedVector<BaseColumns>&& base_structure, const std::vector<RunData>& run_data, const ulint domain) {
        auto run_cols_widths = get_run_cols_widths(run_data);
        if constexpr (IntegratedMoveStructure) {
            auto base_widths = base_structure.get_widths();
            auto widths = get_widths(base_widths, run_cols_widths);

            PackedVector<Columns> final_structure(base_structure.size(), widths);
            for (size_t i = 0; i < final_structure.size(); ++i) {
                auto base_row = base_structure.template get_row(i);
                auto run_row = run_data[i];
                auto row = get_row(base_row, run_row);
                final_structure.set_row(i, row);
            }
            move_structure = MoveStructurePerm(std::move(final_structure), domain);
        } else {
            move_structure = MoveStructurePerm(std::move(base_structure), domain);
            fill_seperated_data(run_data, run_cols_widths);
        }
    }
};

// A wrapper around RunPerm without any run data, essentially just a MoveStructure
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, typename BaseColumns = MoveCols, template<typename, template<typename> class> class MoveStructureType = MoveStructure, template<typename> class TableType = MoveVector>
class MovePermImpl {
protected:
    using RunPermType = RunPermImpl<EmptyRunCols, false, StoreAbsolutePositions, BaseColumns, MoveStructureType, TableType>;
    RunPermType run_perm;
    
public:
    using Position = typename RunPermType::Position;
    
    MovePermImpl() = default;
    
    // Constructor from permutation vector
    MovePermImpl(std::vector<ulint>& permutation, SplitParams split_params = SplitParams()) {
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
        std::vector<std::array<ulint, 0>> empty_run_data(lengths.size());
        run_perm = RunPermType(lengths, interval_permutation, permutation.size(), split_params, empty_run_data);
    }
    
    // Constructor from lengths and interval permutation
    MovePermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, SplitParams split_params = SplitParams()) {
        std::vector<std::array<ulint, 0>> empty_run_data(lengths.size());
        run_perm = RunPermType(lengths, interval_permutation, domain, split_params, empty_run_data);
    }
    
    // Delegate all RunPerm methods
    void first() { run_perm.first(); }
    void last() { run_perm.last(); }
    Position get_position() const { return run_perm.get_position(); }
    void next() { run_perm.next(); }
    void next(ulint steps) { run_perm.next(steps); }
    bool up() { return run_perm.up(); }
    bool down() { return run_perm.down(); }
    ulint get_length() const { return run_perm.get_length(); }
    
    ulint size() const { return run_perm.size(); }
    ulint move_runs() const { return run_perm.move_runs(); }
    ulint permutation_runs() const { return run_perm.permutation_runs(); }
    
    size_t serialize(std::ostream& out) { return run_perm.serialize(out); }
    void load(std::istream& in) { run_perm.load(in); }
};

#endif /* end of include guard: _INTERNAL_RUNPERM_HPP */