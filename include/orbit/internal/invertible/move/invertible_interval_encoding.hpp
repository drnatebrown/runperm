#ifndef _INVERTIBLE_INTERVAL_ENCODING_HPP
#define _INVERTIBLE_INTERVAL_ENCODING_HPP

#include "orbit/common.hpp"
#include "orbit/internal/move/interval_encoding_impl.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include "orbit/internal/ds/packed_vector_aligned.hpp"

namespace orbit {

template<typename int_vector_t = int_vector_aligned>
class invertible_interval_encoding_impl : public interval_encoding_impl<int_vector_t> {
    using base = interval_encoding_impl<int_vector_t>;
public:

private:
    int_vector_t is_fwd_interval;
    int_vector_t is_inv_interval;
};

} // namespace orbit

#endif /* end of include guard: _INVERTIBLE_INTERVAL_ENCODING_HPP */