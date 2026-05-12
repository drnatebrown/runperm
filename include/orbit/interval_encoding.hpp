#ifndef _PUBLIC_INTERVAL_ENCODING_HPP
#define _PUBLIC_INTERVAL_ENCODING_HPP

#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/rlbwt/specializations/rlbwt_interval_encoding.hpp"

namespace orbit {

using interval_encoding = interval_encoding_impl<false, int_vector_aligned>;
using invertible_interval_encoding = interval_encoding_impl<true, int_vector_aligned>;

} // namespace orbit

#endif /* end of include guard: _PUBLIC_INTERVAL_ENCODING_HPP */