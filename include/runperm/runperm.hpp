#ifndef _RUNPERM_HPP
#define _RUNPERM_HPP

#include "common.hpp"
#include "move/move_structure.hpp"
#include "runperm/run_columns.hpp"

constexpr bool DEFAULT_INTEGRATED_MOVE_STRUCTURE = true;
constexpr bool DEFAULT_SPLIT_INTERVALS = true;
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
         bool SplitIntervals = DEFAULT_SPLIT_INTERVALS, // Whether to split intervals if it minimizes table size
         bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, // Whether to store absolute positions instead of interval/offset
         typename BaseColumns = MoveCols,
         template<typename> class TableType = MoveVector,
         template<typename> class PackedType = PackedVector>
class RunPerm : SeperatedDataHolder<RunColsType, IntegratedMoveStructure> {
private:
    static constexpr size_t NumRunCols = static_cast<size_t>(RunColsType::COUNT);
    using SwitchedColumns = SwitchColumns<BaseColumns, StoreAbsolutePositions>;
    // Use RunDataColumns if integrating user data alongside the move structure, otherwise just use the switched columns
    // Use the "fake" enum from RunDataColumns to get the correct columns type for integrated move structure
    using ColumnsType = std::conditional_t<IntegratedMoveStructure, 
        typename RunColsWrapper<RunColsType, SwitchedColumns>::E, 
        SwitchedColumns>;
    MOVE_CLASS_TRAITS(ColumnsType)
    using Table = TableType<ColumnsType>;
    using MoveStructureType = MoveStructure<Table>;

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

    RunPerm(std::vector<ulint> &permutation, const std::vector<std::array<ulint, NumRunCols>> &run_data) {
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
        RunPerm(lengths, interval_permutation, run_data);
    }

    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const std::vector<std::array<ulint, NumRunCols>> &run_data) {
        position = Position(); // should start at 0
        auto run_cols_widths = get_run_cols_widths(run_data);
        PermutationStats stats(lengths);
        
        if constexpr (IntegratedMoveStructure) {
            auto base_widths = MoveStructureType::get_move_widths(stats);
            auto widths = base_widths;
            set_run_cols_widths(widths, run_cols_widths, std::make_index_sequence<NumRunCols>{});
            move_structure = MoveStructureType::construct_with_extended_columns(lengths, interval_permutation, stats, widths, [&] (auto &structure, size_t interval_idx, size_t tbl_idx) {
                set_run_cols_fields(structure, tbl_idx, run_data[interval_idx],
                        std::make_index_sequence<NumRunCols>{});
            });
        } else {
            move_structure = MoveStructureType(lengths, interval_permutation, stats);
            this->run_cols_data = PackedVector<RunCols>(run_data.size(), run_cols_widths);
            for (size_t i = 0; i < run_data.size(); ++i) {
                this->run_cols_data.template set_row(i, run_data[i]);
            }
        }
        // std::optional<ulint> optimal_length_max = std::nullopt;
        // constrexpr if (SplitIntervals) {
        //     optimal_length_max = find_optimal_length_max(lengths, interval_permutation, optimal_length_max);
        // }
    }
            
    // TODO METADATA TO RECTIFY SPLITS?
    // template<auto RunDataField>
    // void set_data(size_t run_id, ulint value) {
    //     storage.template set_data<UserField>(run_id, value);
    // }
    
    // template<auto RunDataField>
    // ulint get_data(size_t run_id) const {
    //     return storage.template get_data<UserField>(run_id);
    // }
    
    void next() { position = move_structure.move(position); }
    // PREVIOUS IF INVERSE SUPPORTED

    void first() { position = move_structure.first(); }
    void last() { position = move_structure.last(); }

    Position get_position() const { return position; }

    ulint size() const { return move_structure.runs(); }
    ulint intervals() const { return orig_intervals; }
    ulint domain() const { return move_structure.size(); }

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
        if constexpr (IntegratedMoveStructure) {
            written_bytes += this->run_cols_data.serialize(os);
        }
        return written_bytes;
    }

    void load(std::istream& is) {
        move_structure.load(is);
        if constexpr (IntegratedMoveStructure) {
            this->run_cols_data.load(is);
        }
    }

private:
    MoveStructureType move_structure;
    Position position;
    size_t orig_intervals; // before splitting, underlying move structure may be larger

    template<typename S, size_t... I>
    static void set_run_cols_fields(S& structure, size_t tbl_idx,
                                const std::array<ulint, NumRunCols>& row,
                                std::index_sequence<I...>) {
        (structure.template set<
            ColsTraits::template run_column<static_cast<RunCols>(I)>()
        >(tbl_idx, row[I]), ...);
    }

    std::array<uchar, NumRunCols> get_run_cols_widths(const std::vector<std::array<ulint, NumRunCols>> &run_data) {
        std::array<ulint, NumRunCols> max_value = {0};
        std::array<uchar, NumRunCols> run_cols_widths = {0};
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

    template<size_t... I>
    static void set_run_cols_widths(std::array<uchar, ColsTraits::NUM_COLS>& widths,
                            const std::array<uchar, NumRunCols>& fw,
                            std::index_sequence<I...>) {
        ((widths[static_cast<size_t>(ColsTraits::template run_column<static_cast<RunCols>(I)>())] = fw[I]), ...);
    }

    // Find the optimal length max for the intervals
    // std::optional<ulint> find_optimal_length_max(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation) {
    //     if constexpr (is_move_table_type()) {
    //         return MAX_VAL(Table::Row::OFFSET_BITS);
    //     } else {
    //         std::array<
    //     }
    // }
};

// A wrapper around RunPerm without any run data, essentially just a MoveStructure
template<bool SplitIntervals = DEFAULT_SPLIT_INTERVALS, bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, 
    typename BaseColumns = MoveCols, template<typename> class TableType = MoveVector, template<typename> class PackedType = PackedVector>
class MovePerm {
private:
    using RunPermType = RunPerm<void, false, SplitIntervals, StoreAbsolutePositions, BaseColumns, TableType, PackedType>;
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
    Position first() const { return run_perm.first(); }
    Position last() const { return run_perm.last(); }
    Position move(Position pos) const { return run_perm.move(pos); }
    Position get_position() const { return run_perm.get_position(); }
    void next() { run_perm.next(); }
    
    ulint size() const { return run_perm.size(); }
    ulint intervals() const { return run_perm.intervals(); }
    ulint domain() const { return run_perm.domain(); }
    
    size_t serialize(std::ostream& out) { return run_perm.serialize(out); }
    void load(std::istream& in) { run_perm.load(in); }
};

#endif /* end of include guard: _RUNPERM_HPP */