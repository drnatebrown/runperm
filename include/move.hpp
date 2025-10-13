// Convenience Header for MoveStructure
#ifndef _MOVE_HPP
#define _MOVE_HPP

#include "internal/common.hpp"
#include "internal/move/move_structure.hpp"

// No simplified interface for MoveStructure, see documentation in include/internal/move/move_structure.hpp

using MoveStructureTbl = MoveStructure<MoveCols, MoveTable>;
using MoveStructureTblIdx = MoveStructure<MoveColsIdx, MoveTable>;
using MoveStructureVec = MoveStructure<MoveCols, MoveVector>;
using MoveStructureVecIdx = MoveStructure<MoveColsIdx, MoveVector>;

#endif /* end of include guard: _MOVE_HPP */