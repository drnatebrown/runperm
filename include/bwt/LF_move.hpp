// #ifndef _LF_move_HH
// #define _LF_move_HH

// #include <fstream>
// #include <iostream>
// #include <optional>
// #include <string>
// #include <vector>
// #include <assert.h>
// #include <cstdio>
// #include "common.hpp"
// #include "ds/move_row.hpp"
// #include "ds/packed_row.hpp"

// using namespace std;

// struct position {
//     ulint interval = 0;
//     ulint offset = 0;
// };

// template <typename row_t = MoveRowStruct>
// class LF_move
// {
// public:
//     LF_move() {}

//     // Expects two files: one with characters (BYTE) and one with lengths (5 BYTES)
//     LF_move(std::ifstream &heads, std::ifstream &lengths)
//     {
//         heads.clear();
//         heads.seekg(0);
//         lengths.clear();
//         lengths.seekg(0);
        
//         LF_runs = vector<row_t>();
//         vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
//         int c_in;
//         ulint i = 0;
//         r = 0;
//         n = 0;
//         while ((c_in = heads.get()) != EOF)
//         {
//             uchar c = static_cast<uchar>(c_in);
//             size_t length = 0;
//             lengths.read((char *)&length, RW_BYTES);
//             if (c <= TERMINATOR) c = TERMINATOR;

//             LF_runs.push_back({c, length, 0, 0});
//             L_block_indices[c].push_back(i++);
//             n+=length;
//         }
//         r = LF_runs.size();

//         compute_table(L_block_indices);
//     }

//     // Expects a file with the BWT (BYTE per character)
//     LF_move(std::ifstream &bwt)
//     {
//         bwt.clear();
//         bwt.seekg(0);
        
//         LF_runs = vector<row_t>();
//         vector<vector<size_t>> L_block_indices = vector<vector<size_t>>(ALPHABET_SIZE);
        
//         uchar last_c;
//         int c_in;
//         ulint run = 0;
//         ulint idx = 0;
//         r = 0;
//         n = 0;
//         size_t length = 0;
//         while ((c_in = bwt.get()) != EOF)
//         {
//             uchar c = static_cast<uchar>(c_in);
//             if (c <= TERMINATOR) c = TERMINATOR;
            
//             if (idx != 0 && c != last_c)
//             {
//                 LF_runs.push_back({last_c, length, 0, 0});
//                 L_block_indices[last_c].push_back(run++); 
//                 n += length;
//                 length = 0;
//             }
//             ++length;
//             ++idx;
//             last_c = c;
//         }
//         // Step for final character
//         LF_runs.push_back({last_c, length, 0, 0});
//         L_block_indices[last_c].push_back(run++);  
//         n+=length;

//         r = LF_runs.size();

//         compute_table(L_block_indices);
//     }

//     const row_t get(size_t i)
//     {
//         assert(i < LF_runs.size());
//         return LF_runs[i];
//     }

//     const row_t get(position p)
//     {
//         assert(p.interval < LF_runs.size());
//         return LF_runs[p.interval];
//     }

//     uchar get_char(size_t i)
//     {
//         return get(i).get_character();
//     }

//     uchar get_char(position p)
//     {
//         return get(p).get_character();
//     }

//     ulint size()
//     {
//         return n;
//     }

//     ulint runs()
//     {
//         return r;
//     }

//     void invert(std::string outfile) 
//     {
//         std::ofstream out(outfile);

//         position pos = {0, 0};

//         uchar c;
//         while((c = get_char(pos)) > TERMINATOR) 
//         {
//             out << c;
//             pos = LF(pos);
//         }
//         out.close();
//     }

//     /*
//      * \param Position (run/offset pair)
//      * \return position (run/offset) of preceding character in text
//      */
//     position LF(position pos)
//     {
//         position next_pos = {get(pos).get_pointer(), get(pos).get_offset() + pos.offset};

//         while (next_pos.offset >= get(next_pos.interval).get_length()) 
//         {
//             next_pos.offset -= get(next_pos.interval++).get_length();
//         }

//         return next_pos;
//     }

//     /* Returns row of largest idx before or at position run with character c */
//     std::optional<position> pred_char(position pos, uchar c)
//     {
//         while (get_char(pos) != c) 
//         {
//             if (pos.run == 0) return std::nullopt;
//             --pos.interval;
//         }
//         pos.offset = get(pos.interval).get_length() - 1;

//         return pos;
//     }

//     /* Returns row of smallest idx after or at position run with character c */
//     std::optional<position> succ_char(position pos, uchar c)
//     {
//         while (get_char(pos) != c)  
//         {
//             if (pos.interval == r - 1) return std::nullopt;
//             ++pos.interval;
//         }
//         pos.offset = 0;

//         return pos;
//     }

//     std::string get_file_extension() const
//     {
//         return ".LF_move";
//     }

//     void bwt_stats()
//     {
//         cout << "Number of BWT equal-letter runs: r = " << r << std::endl;
//         cout << "Length of complete BWT: n = " << n << std::endl;
//         cout << "Rate n/r = " << double(n) / r << std::endl;
//     }

//     /* serialize to the ostream
//     * \param out     the ostream
//     */
//     size_t serialize(std::ostream &out)
//     {
//         size_t written_bytes = 0;

//         out.write((char *)&n, sizeof(n));
//         written_bytes += sizeof(n);

//         out.write((char *)&r, sizeof(r));
//         written_bytes += sizeof(r);

//         const char* data = reinterpret_cast<const char*>(LF_runs.data());
//         size_t size = LF_runs.size() * sizeof(row_t);
//         out.write(data, size);
//         written_bytes += size;

//         return written_bytes;
//     }

//     /* load from the istream
//     * \param in the istream
//     */
//     void load(std::istream &in)
//     {
//         size_t size;

//         in.read((char *)&n, sizeof(n));
//         in.read((char *)&r, sizeof(r));

//         LF_runs = std::vector<row_t>(r);
//         char* data = reinterpret_cast<char*>(LF_runs.data());
//         size = LF_runs.size() * sizeof(row_t);
//         in.read(data, size);
//     }

// protected:
//     ulint n; // Length of BWT
//     ulint r; // Runs of BWT

//     vector<row_t> LF_runs;

//     void compute_table(vector<vector<ulint>> L_block_indices) {\
//         ulint curr_L_num = 0;
//         ulint L_seen = 0;
//         ulint F_seen = 0;
//         for(size_t i = 0; i < L_block_indices.size(); ++i) 
//         {
//             for(size_t j = 0; j < L_block_indices[i].size(); ++j) 
//             {
//                 ulint pos = L_block_indices[i][j];

//                 LF_runs[pos].pointer = curr_L_num;
//                 LF_runs[pos].offset = F_seen - L_seen;

//                 F_seen += get(pos).length;

//                 while (curr_L_num < r && F_seen >= L_seen + get(curr_L_num).length) 
//                 {
//                     L_seen += get(curr_L_num).length;
//                     ++curr_L_num;
//                 }
//             }
//         }
//     }
// };

// #endif /* end of include guard: _LF_move_HH */