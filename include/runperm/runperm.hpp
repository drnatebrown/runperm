#ifndef _RUNPERM_HPP
#define _RUNPERM_HPP

#include "common.hpp"
#include "move.hpp"
#include "runperm/run_data.hpp"

// If we're integrating the run data alongside the move structure, we don't need to store it separately
template <typename RunData, bool IntegratedMoveStructure>
struct FieldDataHolder;
template <typename RunData>
struct FieldDataHolder<RunData, false> { [[no_unique_address]] PackedVector<RunData> field_data; };
template <typename RunData>
struct FieldDataHolder<RunData, true> { /* empty */ };

// TODO InversePermutation, which builds both the forward and inverse move structures if needed
template<typename RunDataType, // Fields to be stored alongside the move structure representing a runny permutation
         bool IntegratedMoveStructure = true, // Whether to pack the run data alongside the move structure
         bool SplitIntervals = true, // Whether to split intervals if it minimizes table size
         bool StoreAbsolutePositions = false, // Whether to store absolute positions instead of interval/offset
         typename BaseColumns = MoveCols,
         template<typename> class TableType = MoveVector,
         template<typename> class PackedType = PackedVector>
class RunPerm : FieldDataHolder<RunDataType, IntegratedMoveStructure> {
    using RunData = RunDataType;
private:
    static constexpr size_t NUM_FIELDS = static_cast<size_t>(RunDataType::COUNT);

    // check if we're using MoveTable
    static constexpr bool is_move_table_type() {
        return std::is_same_v<TableType<void>, MoveTable<void>>;
    }
    // At compile time, check if trying to integrate user data to MoveTable, which uses bitpacked structs
    static_assert(!(IntegratedMoveStructure && is_move_table_type()),
       "Cannot integrate user data with MoveTable. Use MoveVector or set IntegratedMoveStructure=false."
       "MoveTable for integrated data requires specialized implementation!");

    // Switch between relative and absolute columns if needed
    using SwitchedColumns = SwitchColumns<BaseColumns, StoreAbsolutePositions>;
    // Use RunDataColumns if integrating user data alongside the move structure, otherwise just use the switched columns
    using ColumnsHelper = RunDataColumns<RunData, SwitchedColumns>;
    using ColumnsType = std::conditional_t<IntegratedMoveStructure, typename ColumnsHelper::E, SwitchedColumns>;
    MOVE_CLASS_TRAITS(ColumnsType)
    using Table = TableType<Columns>;
    using Permutation = MoveStructure<Table>;
    // using Position = typename Permutation::Position;

    template<size_t... I>
    static void set_user_widths(std::array<uchar, ColsTraits::NUM_COLS>& widths,
                            const std::array<uchar, NUM_FIELDS>& fw,
                            std::index_sequence<I...>) {
        ((widths[static_cast<size_t>(ColumnsHelper::template field<static_cast<RunData>(I)>())] = fw[I]), ...);
    }

    template<typename S, size_t... I>
    static void set_user_fields(S& structure, size_t tbl_idx,
                                const std::array<ulint, NUM_FIELDS>& row,
                                std::index_sequence<I...>) {
        (structure.template set<
            ColumnsHelper::template field<static_cast<RunData>(I)>()
        >(tbl_idx, row[I]), ...);
    }

public:
    using Position = typename Permutation::Position;
    
    RunPerm(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const std::vector<std::array<ulint, NUM_FIELDS>> &run_data) {
        position = Position(); // should start at 0
        auto field_widths = get_field_widths(run_data);
        PermutationStats stats(lengths);
        
        if constexpr (IntegratedMoveStructure) {
            auto base_widths = Permutation::get_move_widths(stats);
            auto widths = base_widths;
            set_user_widths(widths, field_widths, std::make_index_sequence<NUM_FIELDS>{});
            perm = Permutation::construct_with_extended_columns(lengths, interval_permutation, stats, widths, [&] (auto &structure, size_t interval_idx, size_t tbl_idx) {
                for (size_t j = 0; j < NUM_FIELDS; ++j) {
                    set_user_fields(structure, tbl_idx, run_data[interval_idx],
                        std::make_index_sequence<NUM_FIELDS>{});
                }
            });
        } else {
            perm = Permutation(lengths, interval_permutation, stats);
            this->field_data = PackedVector<RunData>(run_data.size(), field_widths);
            for (size_t i = 0; i < run_data.size(); ++i) {
                this->field_data.template set_row(i, run_data[i]);
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
    
    void next() { position = perm.move(position); }
    // PREVIOUS IF INVERSE SUPPORTED

    void first() { position = perm.first(); }
    void last() { position = perm.last(); }

    Position get_position() const { return position; }

    template<RunData Field>
    ulint get() const {
        if constexpr (IntegratedMoveStructure) {
            return perm.template get<ColumnsHelper::template field<Field>()>(position.interval);
        } else {
            return this->field_data.template get<Field>(position.interval);
        }
    }

private:
    Permutation perm;
    Position position;

    std::array<uchar, NUM_FIELDS> get_field_widths(const std::vector<std::array<ulint, NUM_FIELDS>> &run_data) {
        std::array<ulint, NUM_FIELDS> max_value = {0};
        std::array<uchar, NUM_FIELDS> field_widths = {0};
        for (size_t i = 0; i < NUM_FIELDS; ++i) {
            for (size_t j = 0; j < run_data.size(); ++j) {
                max_value[i] = std::max(max_value[i], run_data[j][i]);
            }
        }
        for (size_t i = 0; i < NUM_FIELDS; ++i) {
            field_widths[i] = bit_width(max_value[i]);
        }
        return field_widths;
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
#endif /* end of include guard: _RUNPERM_HPP */