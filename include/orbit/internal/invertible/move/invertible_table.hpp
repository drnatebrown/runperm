#ifndef ORBIT_INVERTIBLE_TABLE_HPP
#define ORBIT_INVERTIBLE_TABLE_HPP

#include "orbit/common.hpp"

namespace orbit {

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

} // namespace orbit

#endif /* end of include guard: ORBIT_INVERTIBLE_TABLE_HPP */