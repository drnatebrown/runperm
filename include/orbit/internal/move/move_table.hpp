#ifndef _MOVE_TABLE_HH
#define _MOVE_TABLE_HH

#include "orbit/common.hpp"
#include "orbit/internal/move/move_row.hpp"
#include "orbit/internal/move/move_columns.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_row.hpp"
#include "orbit/internal/ds/packed_vector.hpp"

#include <cassert>
#include <numeric>

namespace orbit {

template<typename derived, typename columns_t>
struct move_table_interface {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(columns_t)

    template <typename C = columns>
    void set_primary(size_t i, ulint start, ulint length) {
        if constexpr (cols_traits_for<C>::RELATIVE) {
            static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::PRIMARY)>(i, length);
        } else {
            static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::PRIMARY)>(i, start);
        }
    }

    template <typename C = columns>
    std::enable_if_t<cols_traits_for<C>::RELATIVE, void>
    set_length(size_t i, ulint l) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::LENGTH)>(i, l);
    }
    
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, void>
    set_start(size_t i, ulint s) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::START)>(i, s);
    }
    
    void set_pointer(size_t i, ulint p) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::POINTER)>(i, p);
    }
    
    void set_offset(size_t i, ulint o) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::OFFSET)>(i, o);
    }
    
    ulint get_primary(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::PRIMARY)>(i);
    }

    template <typename C = columns>
    std::enable_if_t<cols_traits_for<C>::RELATIVE, ulint>
    get_length(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits_for<C>::LENGTH)>(i);
    }
    
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, ulint>
    get_start(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits_for<C>::START)>(i);
    }   
    
    ulint get_pointer(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::POINTER)>(i);
    }
    
    ulint get_offset(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::OFFSET)>(i);
    }
};

// ============================================= INVERTIBLE =============================================

template<typename derived, typename columns_t>
struct invertible_table_interface {
    // Sets NumCols, Columns, and ColsTraits
    MOVE_CLASS_TRAITS(columns_t)

    template <typename C = columns>
    void set_primary(size_t i, ulint start, ulint length) {
        if constexpr (cols_traits_for<C>::RELATIVE) {
            static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::PRIMARY)>(i, length);
        } else {
            static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::PRIMARY)>(i, start);
        }
    }

    template <typename C = columns>
    std::enable_if_t<cols_traits_for<C>::RELATIVE, void>
    set_length(size_t i, ulint l) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::LENGTH)>(i, l);
    }
    
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, void>
    set_start(size_t i, ulint s) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits_for<C>::START)>(i, s);
    }
    
    void set_pointer_fwd(size_t i, ulint p) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::POINTER_FWD)>(i, p);
    }

    void set_pointer_inv(size_t i, ulint p) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::POINTER_INV)>(i, p);
    }

    void set_fwd_interval(size_t i, bool is_fwd) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::FWD_INTERVAL)>(i, static_cast<ulint>(is_fwd));
    }

    void set_inv_interval(size_t i, bool is_inv) {
        static_cast<derived*>(this)->template set<to_cols(cols_traits::INV_INTERVAL)>(i, static_cast<ulint>(is_inv));
    }
    
    ulint get_primary(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::PRIMARY)>(i);
    }

    template <typename C = columns>
    std::enable_if_t<cols_traits_for<C>::RELATIVE, ulint>
    get_length(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits_for<C>::LENGTH)>(i);
    }
    
    template <typename C = columns>
    std::enable_if_t<!cols_traits_for<C>::RELATIVE, ulint>
    get_start(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits_for<C>::START)>(i);
    }   
    
    ulint get_pointer_fwd(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::POINTER_FWD)>(i);
    }
    
    ulint get_pointer_inv(size_t i) const {
        return static_cast<const derived*>(this)->template get<to_cols(cols_traits::POINTER_INV)>(i);
    }

    bool get_fwd_interval(size_t i) const {
        return static_cast<bool>(static_cast<const derived*>(this)->template get<to_cols(cols_traits::FWD_INTERVAL)>(i));
    }

    bool get_inv_interval(size_t i) const {
        return static_cast<bool>(static_cast<const derived*>(this)->template get<to_cols(cols_traits::INV_INTERVAL)>(i));
    }
};

// ============================================= MOVE TABLE =============================================

template<typename columns_t = move_columns, template<typename, typename> class interface = move_table_interface, typename row_t = typename table_row_for<columns_t>::type>
struct move_table_impl : public interface<move_table_impl<columns_t, interface, row_t>, columns_t> {
    using row = row_t;
    using row_traits = typename row::row_traits;
    // Sets num_cols, columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)
    
    std::vector<row> table;

    move_table_impl() = default;
    move_table_impl(packed_vector<columns> &&vec) {
        row::assert_widths(vec.get_widths());
        
        table = std::vector<row>(vec.size());
        for (size_t i = 0; i < vec.size(); ++i) {
            set_row(i, vec.get_row(i));
        }
    }

    const std::array<uchar, num_cols>& get_widths() const {
        return row::get_widths();
    }

    size_t size() const { return table.size(); }

    template <columns col>
    void set(size_t i, ulint val) {
        table[i].template set<col>(val);
    }

    template <columns col>
    ulint get(size_t i) const {
        return table[i].template get<col>();
    }

    void set_row(size_t i, const std::array<ulint, num_cols>& values) {
        table[i].set(values);
    }

    std::array<ulint, num_cols> get_row(size_t i) const {
        return table[i].get();
    }

    size_t serialize(std::ostream &out) {
        size_t written_bytes = 0;

        size_t tbl_size = table.size();
        out.write((char *)&tbl_size, sizeof(tbl_size));
        written_bytes += sizeof(tbl_size);

        char* data = reinterpret_cast<char*>(table.data());
        size_t size = tbl_size * sizeof(row);
        out.write(data, size);
        written_bytes += size;

        return written_bytes;
    }

    void load(std::istream &in)
    {
        size_t size;
        in.read((char *)&size, sizeof(size));

        table = std::vector<row>(size);
        char* data = reinterpret_cast<char*>(table.data());
        size_t bytes = size * sizeof(row);
        in.read(data, bytes);
    }

    // Widths don't help, the struct is already defined
    static size_t bits_needed(size_t num_rows, std::array<uchar, num_cols> widths) {
        return bytes_to_bits(sizeof(row)) * num_rows;
    }
};

template <typename columns_t = move_columns, template<typename, typename> class interface = move_table_interface>
struct move_vector_impl : public interface<move_vector_impl<columns_t, interface>, columns_t> {
    // Sets num_cols, columns, and cols_traits
    MOVE_CLASS_TRAITS(columns_t)
    
    packed_vector<columns> vec;

    move_vector_impl() = default;
    move_vector_impl(packed_vector<columns> &&vec) : vec(std::move(vec)) {}

    size_t size() const { return vec.size(); }

    template <columns col>
    void set(size_t i, ulint val) {
        vec.template set<col>(i, val);
    }

    template <columns col>
    ulint get(size_t i) const {
        return vec.template get<col>(i);
    }

    void set_row(size_t i, std::array<ulint, num_cols> values) {
        vec.set_row(i, values);
    }

    std::array<ulint, num_cols> get_row(size_t i) const {
        return vec.get_row(i);
    }

    const std::array<uchar, num_cols>& get_widths() const {
        return vec.get_widths();
    }

    size_t serialize(std::ostream &out) {
        return vec.serialize(out);
    }

    void load(std::istream &in)
    {
        vec.load(in);
    }

    // Easy, just sum the widths
    static size_t bits_needed(size_t num_rows, std::array<uchar, num_cols> widths) {
        size_t total_width = std::accumulate(widths.begin(), widths.end(), 0, [](size_t sum, uchar width) {
            return sum + width;
        });
        return total_width * num_rows;
    }
};

template<typename columns_t>
using move_vector = std::conditional_t<resolve_cols_traits<columns_t>::type::INVERTIBLE, 
    move_vector_impl<columns_t, invertible_table_interface>,
    move_vector_impl<columns_t, move_table_interface>>;

template<typename columns_t>
using move_table = std::conditional_t<resolve_cols_traits<columns_t>::type::INVERTIBLE, 
    move_table_impl<columns_t, invertible_table_interface>,
    move_table_impl<columns_t, move_table_interface>>;

} // namespace orbit

#endif // _MOVE_TABLE_HH