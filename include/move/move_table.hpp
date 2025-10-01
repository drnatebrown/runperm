#ifndef _MOVE_TABLE_HH
#define _MOVE_TABLE_HH

#include "common.hpp"
#include "move/move_row.hpp"
#include "move/move_columns.hpp"
#include "ds/packed_vector.hpp"

template<typename Derived, typename ColumnsType>
struct MoveTableInterface {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(ColumnsType)

    template <typename C = Columns>
    void set_primary(size_t i, ulint start, ulint length) {
        if constexpr (MoveColsTraits<C>::RELATIVE) {
            static_cast<Derived*>(this)->template set<to_cols(MoveColsTraits<C>::PRIMARY)>(i, length);
        } else {
            static_cast<Derived*>(this)->template set<to_cols(MoveColsTraits<C>::PRIMARY)>(i, start);
        }
    }

    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::HAS_LENGTH, void>
    set_length(size_t i, ulint l) {
        static_cast<Derived*>(this)->template set<to_cols(MoveColsTraits<C>::LENGTH)>(i, l);
    }
    
    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::HAS_START, void>
    set_start(size_t i, ulint s) {
        static_cast<Derived*>(this)->template set<to_cols(MoveColsTraits<C>::START)>(i, s);
    }
    
    void set_pointer(size_t i, ulint p) {
        static_cast<Derived*>(this)->template set<to_cols(ColsTraits::POINTER)>(i, p);
    }
    
    void set_offset(size_t i, ulint o) {
        static_cast<Derived*>(this)->template set<to_cols(ColsTraits::OFFSET)>(i, o);
    }
    
    ulint get_primary(size_t i) const {
        return static_cast<const Derived*>(this)->template get<to_cols(ColsTraits::PRIMARY)>(i);
    }

    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::HAS_LENGTH, ulint>
    get_length(size_t i) const {
        return static_cast<const Derived*>(this)->template get<to_cols(MoveColsTraits<C>::LENGTH)>(i);
    }
    
    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::HAS_START, ulint>
    get_start(size_t i) const {
        return static_cast<const Derived*>(this)->template get<to_cols(MoveColsTraits<C>::START)>(i);
    }
    
    ulint get_pointer(size_t i) const {
        return static_cast<const Derived*>(this)->template get<to_cols(ColsTraits::POINTER)>(i);
    }
    
    ulint get_offset(size_t i) const {
        return static_cast<const Derived*>(this)->template get<to_cols(ColsTraits::OFFSET)>(i);
    }
};

template<typename Row = MoveRow<>>
struct MoveTable : public MoveTableInterface<MoveTable<Row>, typename Row::Columns> {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(typename Row::Columns)
    
    std::vector<Row> table;

    MoveTable() = default;
    MoveTable(PackedVector<Columns> &vec) {
        Row::assert_widths(vec.get_widths());
        
        table = std::vector<Row>(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            set_row(i, vec.get_row(i));
        }
    }

    size_t size() const { return table.size(); }

    template <Columns Col>
    void set(size_t i, ulint val) {
        table[i].set<Col>(val);
    }

    template <Columns Col>
    ulint get(size_t i) const {
        return table[i].template get<Col>();
    }

    void set_row(size_t i, const std::array<ulint, NumCols>& values) {
        table[i].set(values);
    }

    std::array<ulint, NumCols> get_row(size_t i) const {
        return table[i].get();
    }

    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        size_t tbl_size = table.size();
        out.write((char *)&tbl_size, sizeof(tbl_size));
        written_bytes += sizeof(tbl_size);

        char* data = reinterpret_cast<char*>(table.data());
        size_t size = tbl_size * sizeof(Row);
        out.write(data, size);
        written_bytes += size;

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t size;
        in.read((char *)&size, sizeof(size));

        table = std::vector<Row>(size);
        char* data = reinterpret_cast<char*>(table.data());
        size_t bytes = size * sizeof(Row);
        in.read(data, bytes);
    }

    // Widths don't help, the struct is already defined
    static size_t bits_needed(size_t num_rows, std::array<uchar, NumCols> widths) {
        return BYTES_TO_BITS(sizeof(Row)) * num_rows;
    }
};

template <typename ColumnsType = MoveCols>
struct MoveVector : public MoveTableInterface<MoveVector<ColumnsType>, ColumnsType> {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(ColumnsType)
    
    PackedVector<Columns> vec;

    MoveVector() = default;
    MoveVector(PackedVector<Columns> &vec) {
        this->vec = vec;
    }

    size_t size() const { return vec.size(); }

    template <Columns Col>
    void set(size_t i, ulint val) {
        vec.set<Col>(i, val);
    }

    template <Columns Col>
    ulint get(size_t i) const {
        return vec.template get<Col>(i);
    }

    void set_row(size_t i, std::array<ulint, NumCols> values) {
        vec.set_row(i, values);
    }

    std::array<ulint, NumCols> get_row(size_t i) const {
        return vec.get_row(i);
    }

    size_t serialize(std::ostream &out) {
        return vec.serialize(out);
    }

    void load(std::istream &in)
    {
        vec.load(in);
    }

    // Easy, just sum the widths
    static size_t bits_needed(size_t num_rows, std::array<uchar, NumCols> widths) {
        size_t total_width = std::accumulate(widths.begin(), widths.end(), 0, [](size_t sum, uchar width) {
            return sum + width;
        });
        return total_width * num_rows;
    }
};

template<typename Columns>
using MoveTableFor = MoveTable<MoveRow<Columns>>;
template<typename Columns>
using MoveVectorFor = MoveVector<Columns>;

using MoveTableIdx = MoveTableFor<MoveColsIdx>;
using MoveVectorIdx = MoveVectorFor<MoveColsIdx>;

#endif // _MOVE_TABLE_HH