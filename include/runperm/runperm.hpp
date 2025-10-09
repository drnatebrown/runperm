#ifndef _RUNPERM_HPP
#define _RUNPERM_HPP

#include "common.hpp"
#include "move/move_structure.hpp"
#include "runperm/run_columns.hpp"

constexpr bool DEFAULT_INTEGRATED_MOVE_STRUCTURE = false;
constexpr bool DEFAULT_STORE_ABSOLUTE_POSITIONS = false;

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
         typename BaseColumns = MoveCols,
         template<typename> class TableType = MoveVector,
         template<typename> class StructureType = MoveStructure,
         template<typename> class PackedType = PackedVector>
class RunPerm : SeperatedDataHolder<RunColsType, IntegratedMoveStructure> {
private:
    static constexpr size_t NumRunCols = static_cast<size_t>(RunColsType::COUNT);
    static constexpr size_t NumBaseCols = static_cast<size_t>(BaseColumns::COUNT);
    // SwitchedColumns is the base colunns using the correct relative/absolute indexing
    using SwitchedColumns = SwitchColumns<BaseColumns, StoreAbsolutePositions>;
    // Use RunDataColumns if integrating user data alongside the move structure, otherwise just use the switched columns
    // Use the "fake" enum from RunDataColumns to get the correct columns type for integrated move structure
    using ColumnsType = std::conditional_t<IntegratedMoveStructure, 
        typename RunColsWrapper<RunColsType, SwitchedColumns>::E, 
        SwitchedColumns>;
    MOVE_CLASS_TRAITS(ColumnsType)
    using Table = TableType<ColumnsType>;
    using MoveStructureType = StructureType<Table>;
    using MoveStructureBase = StructureType<TableType<SwitchedColumns>>;

public:
    using RunCols = RunColsType;
    using Position = typename MoveStructureType::Position;

    // check if we're using MoveTable
    static constexpr bool is_move_table_type() {
        return std::is_same_v<TableType<void>, MoveTable<void>>;
    }

    // At compile time, check if trying to integrate user data to MoveTable, which uses bitpacked structs
    static_assert(!(IntegratedMoveStructure && is_move_table_type()),
    "Cannot integrate user data with MoveTable. Use MoveVector or set IntegratedMoveStructure=false."
    "MoveTable for integrated data requires specialized implementation!");


    RunPerm() = default;

    // By default, just copy the run data for the original interval if the move structure intervals have been split
    /**
     * lengths -> length of each interval which permutes contiguously
     * interval_permutation -> permutation position of the first position of each interval
     * domain -> domain of the permutation, i.e. a permutatation over 1..n has domain n
     * run_data -> run data for each interval, the size of this vector should be the same as the number of intervals
     */
    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, const std::vector<std::array<ulint, NumRunCols>> &run_data)
        : RunPerm(lengths, interval_permutation, domain,
            [&run_data](ulint orig_interval, ulint orig_interval_length, ulint new_offset_from_orig_start, ulint new_length) {
                return run_data[orig_interval];
            }
        ){}

    // Advanced constructor for users who want to specify how to set the run data for each column type by passing a function which is called with:
    // the original interval index, the original interval length, the absolute offset of the first element in the new interval, and the new interval length
    // and returns the values to be stored in the new interval (run data)
    // TODO add skip if no skipping happened
    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, std::function<std::array<ulint, NumRunCols>(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        assert(lengths.size() == interval_permutation.size());
        orig_intervals = lengths.size();
        position = Position(); // should start at 0
        
        PackedVector<SwitchedColumns> base_structure = MoveStructureBase::find_structure(lengths, interval_permutation, domain);
        std::vector<std::array<ulint, NumRunCols>> final_run_data = extend_run_data(lengths, domain, base_structure, get_run_cols_data);
        auto run_cols_widths = get_run_cols_widths(final_run_data);

        if constexpr (IntegratedMoveStructure) {
            auto base_widths = base_structure.get_widths();
            auto widths = get_widths(base_widths, run_cols_widths);

            PackedVector<Columns> final_structure(base_structure.size(), widths);
            for (size_t i = 0; i < final_structure.size(); ++i) {
                auto base_row = base_structure.template get_row(i);
                auto run_row = final_run_data[i];
                auto row = get_row(base_row, run_row);
                final_structure.set_row(i, row);
            }
            move_structure = MoveStructureType(std::move(final_structure), domain);
        } else {
            move_structure = MoveStructureType(std::move(base_structure), domain);
            
            this->run_cols_data = PackedVector<RunCols>(lengths.size(), run_cols_widths);
            for (size_t i = 0; i < final_run_data.size(); ++i) {
                this->run_cols_data.template set_row(i, final_run_data[i]);
            }
        }
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

    size_t serialize(std::ostream& os) {
        size_t written_bytes = 0;
        written_bytes += move_structure.serialize(os);
        if constexpr (!IntegratedMoveStructure) {
            written_bytes += this->run_cols_data.serialize(os);
        }
        return written_bytes;
    }

    void load(std::istream& is) {
        move_structure.load(is);
        if constexpr (!IntegratedMoveStructure) {
            this->run_cols_data.load(is);
        }
    }

private:
    MoveStructureType move_structure;
    Position position;
    size_t orig_intervals; // before splitting, underlying move structure may be larger


    std::vector<std::array<ulint, NumRunCols>> extend_run_data(const std::vector<ulint>& lengths, const ulint domain, const PackedVector<SwitchedColumns>& structure, std::function<std::array<ulint, NumRunCols>(ulint, ulint, ulint, ulint)> get_run_cols_data) {
        std::vector<std::array<ulint, NumRunCols>> final_run_data(structure.size());
        auto get_structure_length = [&](size_t idx) {
            if constexpr (StoreAbsolutePositions) {
                return (idx == structure.size() - 1) 
                    ? domain - structure.template get<ColsTraits::PRIMARY>(idx)
                    : structure.template get<ColsTraits::PRIMARY>(idx + 1) - structure.template get<ColsTraits::PRIMARY>(idx);
            } else {
                return structure.template get<ColsTraits::PRIMARY>(idx);
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

    std::array<uchar, NumRunCols> get_run_cols_widths(const std::vector<std::array<ulint, NumRunCols>> &run_data) {
        std::array<ulint, NumRunCols> max_value = {};
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

    std::array<ulint, NumCols> get_row(const std::array<ulint, NumBaseCols>& base_row, const std::array<ulint, NumRunCols>& run_row) {
        std::array<ulint, NumCols> row;
        for (size_t i = 0; i < NumBaseCols; ++i) {
            row[i] = base_row[i];
        }
        for (size_t i = 0; i < NumRunCols; ++i) {
            row[NumBaseCols + i] = run_row[i];
        }
        return row;
    }
};

// A wrapper around RunPerm without any run data, essentially just a MoveStructure
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, typename BaseColumns = MoveCols, template<typename> class TableType = MoveVector, template<typename> class StructureType = MoveStructure, template<typename> class PackedType = PackedVector>
class MovePerm {
private:
    using RunPermType = RunPerm<EmptyRunCols, false, StoreAbsolutePositions, BaseColumns, TableType, StructureType, PackedType>;
    RunPermType run_perm;
    
public:
    using Position = typename RunPermType::Position;
    
    MovePerm() = default;
    
    // Constructor from permutation vector
    MovePerm(std::vector<ulint>& permutation) {
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
        std::vector<std::array<ulint, 0>> empty_run_data(lengths.size());
        run_perm = RunPermType(lengths, interval_permutation, empty_run_data);
    }
    
    // Constructor from lengths and interval permutation
    MovePerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation) {
        std::vector<std::array<ulint, 0>> empty_run_data(lengths.size());
        run_perm = RunPermType(lengths, interval_permutation, empty_run_data);
    }
    
    // Delegate all RunPerm methods
    void first() { run_perm.first(); }
    void last() { run_perm.last(); }
    void move(Position pos) const { run_perm.move(pos); }
    Position get_position() const { return run_perm.get_position(); }
    void next() { run_perm.next(); }
    
    ulint size() const { return run_perm.size(); }
    ulint move_runs() const { return run_perm.move_runs(); }
    ulint permutation_runs() const { return run_perm.permutation_runs(); }
    
    size_t serialize(std::ostream& out) { return run_perm.serialize(out); }
    void load(std::istream& in) { run_perm.load(in); }
};

#endif /* end of include guard: _RUNPERM_HPP */