#ifndef _INTERNAL_MOVEPERM_HPP
#define _INTERNAL_MOVEPERM_HPP

#include "internal/common.hpp"
#include "internal/move/move_structure.hpp"
#include "internal/runperm/move_columns.hpp"

constexpr bool DEFAULT_STORE_ABSOLUTE_POSITIONS = false;

/* ============================================= Advanced Implemenation ============================================= */
template<bool StoreAbsolutePositions = DEFAULT_STORE_ABSOLUTE_POSITIONS, // Whether to store absolute positions instead of interval/offset
         typename BaseColumnsType = MoveCols,
         template<typename, template<typename> class> class MoveStructureType = MoveStructure,
         template<typename> class TableType = MoveVector>
         // TODO need PackedType option?
class MovePermImpl {
protected:
    static constexpr size_t NumBaseCols = static_cast<size_t>(BaseColumnsType::COUNT);
    using BaseColumns = SwitchColumns<BaseColumnsType, StoreAbsolutePositions>;
    using MoveStructureBase = MoveStructureType<BaseColumns, TableType>;
    MOVE_CLASS_TRAITS(BaseColumns)
public:
    using Position = typename MoveStructureBase::Position;
    
    MovePermImpl() = default;
    
    // Constructor from permutation vector
    MovePermImpl(std::vector<ulint>& permutation, SplitParams split_params = SplitParams()) {
        auto [lengths, interval_permutation] = get_permutation_intervals(permutation);
        assert(lengths.size() == interval_permutation.size());
        orig_intervals = lengths.size();
        move_structure = MoveStructureBase::find_structure(lengths, interval_permutation, permutation.size(), split_params);
    }

    MovePermImpl(const std::vector<ulint>& lengths, const std::vector<ulint>& interval_permutation, const ulint domain, SplitParams split_params = SplitParams()) {
        assert(lengths.size() == interval_permutation.size());
        orig_intervals = lengths.size();
        move_structure = MoveStructureBase::find_structure(lengths, interval_permutation, domain, split_params);
    }

    Position first() { return move_structure.first(); }
    Position last() { return move_structure.last(); }

    Position next(Position position) { return move_structure.move(position); }
    Position next(Position position, ulint steps) {
        for (ulint i = 0; i < steps; ++i) {
            position = move_structure.move(position);
        }
        return position;
    }

    // Set position to interval above in underlying move structure, or circularly wrap to bottom if already at top
    Position up(Position position) {
        if (position.interval == 0)
        {
            return last();
        }
        --position.interval;
        position.offset = move_structure.get_length(position.interval) - 1;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval) + position.offset;
        }
        return position;
    }

    // Set position to interval below in underlying move structure, or nothing if already at bottom
    Position down(Position position) {
        if (position.interval == move_structure.runs() - 1)
        {
            return first();
        }
        ++position.interval;
        position.offset = 0;
        if constexpr (StoreAbsolutePositions) {
            position.idx = move_structure.get_start(position.interval);
        }
        return position;
    }

    ulint get_length(ulint interval) const {
        return move_structure.get_length(interval);
    }
    ulint get_length(Position position) const {
        return get_length(position.interval);
    }
    
    ulint size() const { return move_structure.size(); }
    ulint move_runs() const { return move_structure.runs(); }
    ulint permutation_runs() const { return orig_intervals; }
    
    size_t serialize(std::ostream& os) {
        size_t written_bytes = 0;

        os.write((char *)&orig_intervals, sizeof(orig_intervals));
        written_bytes += sizeof(orig_intervals);

        written_bytes += move_structure.serialize(os);
        return written_bytes;
    }

    void load(std::istream& is) {
        is.read((char *)&orig_intervals, sizeof(orig_intervals));

        move_structure.load(is);
    }

private:
    MoveStructureBase move_structure;
    size_t orig_intervals; // before splitting, underlying move structure may be larger
};

#endif /* end of include guard: _INTERNAL_MOVEPERM_HPP */