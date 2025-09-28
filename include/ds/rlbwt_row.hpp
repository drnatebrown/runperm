// #include "common.hpp"

// using namespace std;

// // Supports π = LF or π = FL
// enum class RLBWT_Col {
//     CHARACTER, // RLBWT Character
//     LENGTH, // RLBWT Length
//     POINTER, // Where BWT[i] is run head, interval is the run containing π(i)
//     OFFSET // Offset of π(i) into the run that contains it
// };

// // Row of the RLBWT table (bitpacked)
// // Bytes set in common.hpp
// template <
//     size_t CHARACTER_BITS = BYTES_TO_BITS(CHARACTER_BYTES),
//     size_t LENGTH_BITS = BYTES_TO_BITS(LENGTH_BYTES),
//     size_t POINTER_BITS = BYTES_TO_BITS(POINTER_BYTES),
//     size_t OFFSET_BITS = BYTES_TO_BITS(OFFSET_BYTES)>
// struct RLBWT_RowStruct {
//     uchar character : CHARACTER_BITS;
//     ulint length : LENGTH_BITS;
//     ulint pointer : POINTER_BITS;
//     ulint offset : OFFSET_BITS;

//     RLBWT_RowStruct() = default;
//     RLBWT_RowStruct(const std::array<ulint, RLBWT_Col::NUM_COLS>& values)
//         : character(values[RLBWT_Col::CHARACTER]), length(values[RLBWT_Col::LENGTH]), pointer(values[RLBWT_Col::POINTER]), offset(values[RLBWT_Col::OFFSET]) {}
//     RLBWT_RowStruct(uchar c, ulint l, ulint p, ulint o) : character(c), length(l), pointer(p), offset(o) {}

//     template <RLBWT_Col Col>
//     uint64_t get() const {
//         if constexpr (Col == RLBWT_Col::CHARACTER) return character;
//         else if constexpr (Col == RLBWT_Col::LENGTH) return length;
//         else if constexpr (Col == RLBWT_Col::POINTER) return pointer;
//         else if constexpr (Col == RLBWT_Col::OFFSET) return offset;
//     }

//     uint64_t get_character() const { return character; }
//     uint64_t get_length() const { return length; }
//     uint64_t get_pointer() const { return pointer; }
//     uint64_t get_offset() const { return offset; }
// } __attribute__((packed));

// struct RLBWT_RowPacked : PackedRow<RLBWT_Col> {
//     RLBWT_RowPacked() = default;
//     RLBWT_RowPacked(const std::array<ulint, RLBWT_Col::NUM_COLS>& values)
//         : PackedRow<RLBWT_Col>(values) {}
//     RLBWT_RowPacked(uchar c, ulint l, ulint p, ulint o)
//         : PackedRow<RLBWT_Col>(std::array<ulint, RLBWT_Col::NUM_COLS>{c, l, p, o}) {}

//     template <RLBWT_Col Col>
//     uint64_t get() const {
//         return PackedRow<RLBWT_Col>::get<Col>();
//     }

//     uint64_t get_character() const { return get<RLBWT_Col::CHARACTER>(); }
//     uint64_t get_length() const { return get<RLBWT_Col::LENGTH>(); }
//     uint64_t get_pointer() const { return get<RLBWT_Col::POINTER>(); }
//     uint64_t get_offset() const { return get<RLBWT_Col::OFFSET>(); }
// };