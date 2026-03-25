// Convenience Header for move_structure
#ifndef _MOVE_HPP
#define _MOVE_HPP

#include "orbit/common.hpp"
#include "orbit/interval_encoding.hpp"
#include "orbit/internal/move/move_structure_impl.hpp"

namespace orbit {

// No simplified interface for move_structure, see documentation in include/internal/move/move_structure.hpp

using move_structure_tbl = move_structure<move_columns, move_table>;
using move_structure_tbl_idx = move_structure<move_columns_idx, move_table>;
using move_structure_vec = move_structure<move_columns, move_vector>;
using move_structure_vec_idx = move_structure<move_columns_idx, move_vector>;

} // namespace orbit``

#endif /* end of include guard: _MOVE_HPP */