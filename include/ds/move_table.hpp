#ifndef _MOVE_TABLE_HH
#define _MOVE_TABLE_HH

#include "common.hpp"
#include "move_row.hpp"
#include "move_columns.hpp"
#include "packed_vector.hpp"

template<typename Derived, typename Columns>
struct MoveTableInterface {
    using Columns = Columns;

    template <typename C = Columns>
    void set_primary(size_t i, ulint start, ulint length) {
        if constexpr (MoveColsTraits<C>::IS_LENGTH) {
            set_length(i, length);
        } else if constexpr (MoveColsTraits<C>::IS_START) {
            set_start(i, start);
        }
    }

    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_LENGTH, void>
    set_length(size_t i, ulint l) {
        static_cast<Derived*>(this)->template set<MoveColsTraits<C>::PRIMARY>(i, l);
    }
    
    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_START, void>
    set_start(size_t i, ulint s) {
        static_cast<Derived*>(this)->template set<MoveColsTraits<C>::PRIMARY>(i, s);
    }
    
    void set_pointer(size_t i, ulint p) {
        static_cast<Derived*>(this)->template set<Columns::POINTER>(i, p);
    }
    
    void set_offset(size_t i, ulint o) {
        static_cast<Derived*>(this)->template set<Columns::OFFSET>(i, o);
    }
    
    ulint get_primary(size_t i) const {
        return static_cast<const Derived*>(this)->template get<MoveColsTraits<Columns>::PRIMARY>(i);
    }

    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_LENGTH, ulint>
    get_length(size_t i) const {
        return static_cast<const Derived*>(this)->template get<MoveColsTraits<C>::PRIMARY>(i);
    }
    
    template <typename C = Columns>
    std::enable_if_t<MoveColsTraits<C>::IS_START, ulint>
    get_start(size_t i) const {
        return static_cast<const Derived*>(this)->template get<MoveColsTraits<C>::PRIMARY>(i);
    }
    
    ulint get_pointer(size_t i) const {
        return static_cast<const Derived*>(this)->template get<Columns::POINTER>(i);
    }
    
    ulint get_offset(size_t i) const {
        return static_cast<const Derived*>(this)->template get<Columns::OFFSET>(i);
    }
};

template<typename Row = MoveRow<>>
struct MoveTable : public MoveTableInterface<MoveTable<Row>, typename Row::Columns> {
    using Columns = typename Row::Columns;
    std::vector<Row> table;

    MoveTable() = default;
    MoveTable(PackedVector<Columns> &vec) {
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
        return table[i].get<Col>();
    }

    void set_row(size_t i, const std::array<ulint, Columns::NUM_COLS>& values) {
        table[i].set(values);
    }

    std::array<ulint, Columns::NUM_COLS> get_row(size_t i) const {
        return table[i].get();
    }

    size_t serialize(std::ostream &out)
    {
        size_t written_bytes = 0;

        out.write((char *)&table.size(), sizeof(table.size()));
        written_bytes += sizeof(table.size());

        char* data = reinterpret_cast<char*>(table.data());
        size_t size = table.size() * sizeof(Row);
        out.write(data, size);
        written_bytes += size;

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t size;
        in.read((char *)&size, sizeof(size));

        table = std::vector<Row>(size);
        const char* data = reinterpret_cast<const char*>(table.data());
        size_t bytes = size * sizeof(Row);
        in.read(data, bytes);
    }
};

template <typename Columns = MoveCols>
struct MoveVector : public MoveTableInterface<MoveVector<Columns>, Columns> {
    using Columns = Columns;
    PackedVector<Columns> vec;

    MoveVector() = default;
    MoveVector(PackedVector<Columns> &vec) {
        this->vec = std::move(vec);
    }

    size_t size() const { return vec.size(); }

    template <Columns Col>
    void set(size_t i, ulint val) {
        vec.set<Col>(i, val);
    }

    template <Columns Col>
    ulint get(size_t i) const {
        return vec.get<Col>(i);
    }

    void set_row(size_t i, std::array<ulint, Columns::NUM_COLS> values) {
        vec.set_row(i, values);
    }

    std::array<ulint, Columns::NUM_COLS> get_row(size_t i) const {
        return vec.get_row(i);
    }

    size_t serialize(std::ostream &out)
    {
        return vec.serialize(out);
    }

    void load(std::istream &in)
    {
        vec.load(in);
    }
};

template<typename Columns>
using MoveTableFor = MoveTable<MoveRow<Columns>>;
template<typename Columns>
using MoveVectorFor = MoveVector<Columns>;

using MoveTableIdx = MoveTableFor<MoveColsIdx>;
using MoveVectorIdx = MoveVectorFor<MoveColsIdx>;

#endif // _MOVE_TABLE_HH