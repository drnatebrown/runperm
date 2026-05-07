#ifndef _RLBWT_HELPERS_HPP
#define _RLBWT_HELPERS_HPP

#include "orbit/common.hpp"

namespace orbit::rlbwt {

inline std::pair<std::vector<uchar>, std::vector<ulint>> bwt_to_rlbwt(const std::vector<uchar> &bwt_chars) {
    std::vector<uchar> rlbwt_chars;
    std::vector<ulint> rlbwt_run_lengths;
    for (size_t i = 0; i < bwt_chars.size(); ++i) {
        if (i == 0 || bwt_chars[i] != bwt_chars[i - 1]) {
            rlbwt_chars.push_back(bwt_chars[i]);
            rlbwt_run_lengths.push_back(1);
        } else {
            ++rlbwt_run_lengths.back();
        }
    }
    return {rlbwt_chars, rlbwt_run_lengths};
}

template<typename container1_t, typename container2_t>
inline std::tuple<std::vector<size_t>, std::vector<size_t>, size_t, ulint> get_LF_char_counts(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths) {
    assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

    std::vector<size_t> char_count(MAX_ALPHABET_SIZE, 0);
    std::vector<size_t> head_ranks(rlbwt_heads.size(), 0);
    size_t bwt_length = 0;
    ulint max_length = 0;
    for (size_t i = 0; i < rlbwt_heads.size(); i++)
    {
        uchar c = rlbwt_heads[i];
        ulint length = rlbwt_run_lengths[i];
        head_ranks[i] = char_count[c];
        char_count[c] += length;
        bwt_length+=length;
        max_length = std::max(max_length, length);
    }
    return {char_count, head_ranks, bwt_length, max_length};
}

template<typename container1_t, typename container2_t>
inline std::tuple<std::vector<size_t>, size_t, ulint> get_LF_head_counts(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths) {
    assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

    std::vector<size_t> head_counts(MAX_ALPHABET_SIZE, 0);
    ulint max_length = 0;
    size_t bwt_length = 0;
    for (size_t i = 0; i < rlbwt_heads.size(); i++)
    {
        uchar c = rlbwt_heads[i];
        ulint length = rlbwt_run_lengths[i];
        ++head_counts[c];
        bwt_length+=length;
        max_length = std::max(max_length, length);
    }
    return {head_counts, bwt_length, max_length};
}

template<typename containter1_t, typename int_vector_t = int_vector_aligned>
inline int_vector_t get_LF_images(const containter1_t& rlbwt_heads, const std::vector<size_t> &char_count, const std::vector<size_t> &head_ranks, size_t bwt_length) {
    int_vector_t images(head_ranks.size(), bit_width(bwt_length - 1));
    
    std::vector<size_t> C_array(char_count.size(), 0);
    size_t seen = 0;
    for (size_t i = 0; i < char_count.size(); i++) {
        C_array[i] = seen;
        seen += char_count[i];
    }

    for (size_t i = 0; i < rlbwt_heads.size(); i++) {
        uchar c = rlbwt_heads[i];
        images[i] = C_array[c] + head_ranks[i];
    }
    return images;
}

template<typename container1_t, typename int_vector_t = int_vector_aligned>
inline int_vector_t get_LF_img_rank_inv(const container1_t& rlbwt_heads, const std::vector<size_t> &head_counts) {
    int_vector_t img_rank_inv(rlbwt_heads.size(), bit_width(rlbwt_heads.size() - 1));
    std::vector<size_t> seen_heads(head_counts.size(), 0);

    std::vector<size_t> C_head_array(head_counts.size(), 0);
    size_t seen = 0;
    for (size_t i = 0; i < head_counts.size(); i++) {
        C_head_array[i] = seen;
        seen += head_counts[i];
    }

    for (size_t i = 0; i < rlbwt_heads.size(); i++) {
        uchar c = rlbwt_heads[i];
        img_rank_inv[C_head_array[c] + seen_heads[c]] = i;
        ++seen_heads[c];
    }

    return img_rank_inv;
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline int_vector_t rlbwt_to_lf_img_rank_inv(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [head_counts, n, max_length_seen] = get_LF_head_counts(rlbwt_heads, rlbwt_run_lengths);
    auto img_rank_inv = get_LF_img_rank_inv(rlbwt_heads, head_counts);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return img_rank_inv;
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline int_vector_t rlbwt_to_lf_images(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [char_counts, head_ranks, n, max_length_seen] = get_LF_char_counts(rlbwt_heads, rlbwt_run_lengths);
    auto images = get_LF_images(rlbwt_heads, char_counts, head_ranks, n);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return images;
}

template<typename container1_t, typename container2_t>
inline std::tuple<std::vector<size_t>, std::vector<std::vector<std::pair<size_t, size_t>>>, size_t, ulint> get_FL_head_counts(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths) {
    assert(rlbwt_heads.size() == rlbwt_run_lengths.size());

    std::vector<size_t> head_counts(MAX_ALPHABET_SIZE, 0);
    // Holds the lengths of runs in F, and the origin run in BWT (i.e., run containing the FL(i) for i at run heads)
    std::vector<std::vector<std::pair<size_t, size_t>>> F_lens_and_origin_run(MAX_ALPHABET_SIZE);
    size_t bwt_length = 0;
    ulint max_length = 0;
    for (size_t i = 0; i < rlbwt_heads.size(); i++)
    {
        uchar c = rlbwt_heads[i];
        ulint length = rlbwt_run_lengths[i];
        F_lens_and_origin_run[c].push_back({length, i});
        ++head_counts[c];
        bwt_length+=length;
        max_length = std::max(max_length, length);
    }
    return {head_counts, F_lens_and_origin_run, bwt_length, max_length};
}

template<typename int_vector_t = int_vector_aligned>
inline std::tuple<std::vector<uchar>, int_vector_t, int_vector_t> get_FL_runs_and_images(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins, size_t bwt_length, ulint max_length) {
    std::vector<uchar> F_heads(runs);
    int_vector_t F_lens(runs, bit_width(max_length));
    int_vector_t images(runs, bit_width(bwt_length - 1));

    size_t curr_run = 0;
    for (size_t c = 0; c < F_lens_and_origins.size(); ++c) {
        for (size_t j = 0; j < F_lens_and_origins[c].size(); ++j) {
            F_heads[curr_run] = c;
            F_lens[curr_run] = F_lens_and_origins[c][j].first;
            images[curr_run] = F_lens_and_origins[c][j].second;
            curr_run++;
        }
    }
    return {F_heads, F_lens, images};
}

template<typename int_vector_t = int_vector_aligned>
inline std::tuple<std::vector<uchar>, int_vector_t, int_vector_t> get_FL_runs_and_img_rank_inv(const size_t runs, const std::vector<std::vector<std::pair<size_t, size_t>>> &F_lens_and_origins, ulint max_length) {
    std::vector<uchar> F_heads(runs);
    int_vector_t F_lens(runs, bit_width(max_length));
    int_vector_t F_img_rank_inv(runs, bit_width(runs - 1));

    size_t curr_run = 0;
    for (size_t c = 0; c < F_lens_and_origins.size(); ++c) {
        for (size_t j = 0; j < F_lens_and_origins[c].size(); ++j) {
            F_heads[curr_run] = c;
            F_lens[curr_run] = F_lens_and_origins[c][j].first;
            F_img_rank_inv[F_lens_and_origins[c][j].second] = curr_run;
            curr_run++;
        }
    }
    return {F_heads, F_lens, F_img_rank_inv};
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline std::tuple<std::vector<uchar>, int_vector_t, int_vector_t> rlbwt_to_fl_runs_and_img_rank_inv(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [head_counts, F_lens_and_origin_run, n, max_length_seen] = get_FL_head_counts(rlbwt_heads, rlbwt_run_lengths);
    auto [F_heads, F_lens, F_img_rank_inv] = get_FL_runs_and_img_rank_inv(rlbwt_heads.size(), F_lens_and_origin_run, max_length_seen);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return {F_heads, F_lens, F_img_rank_inv};
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline std::pair<int_vector_t, int_vector_t> rlbwt_to_fl_img_rank_inv(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [head_counts, F_lens_and_origin_run, n, max_length_seen] = get_FL_head_counts(rlbwt_heads, rlbwt_run_lengths);
    auto [_, F_lens, F_img_rank_inv] = get_FL_runs_and_img_rank_inv(rlbwt_heads.size(), F_lens_and_origin_run, max_length_seen);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return {F_lens, F_img_rank_inv};
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline std::tuple<std::vector<uchar>, int_vector_t, int_vector_t> rlbwt_to_fl_runs_and_images(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [head_counts, F_lens_and_origin_run, n, max_length_seen] = get_FL_head_counts(rlbwt_heads, rlbwt_run_lengths);
    auto [F_heads, F_lens, F_images] = get_FL_runs_and_images(rlbwt_heads.size(), F_lens_and_origin_run, n, max_length_seen);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return {F_heads, F_lens, F_images};
}

template<typename container1_t, typename container2_t, typename int_vector_t = int_vector_aligned>
inline std::pair<int_vector_t, int_vector_t> rlbwt_to_fl_images(const container1_t& rlbwt_heads, const container2_t& rlbwt_run_lengths, size_t* domain = nullptr, ulint* max_length = nullptr) {
    auto [head_counts, F_lens_and_origin_run, n, max_length_seen] = get_FL_head_counts(rlbwt_heads, rlbwt_run_lengths);
    auto [_, F_lens, F_images] = get_FL_runs_and_images(rlbwt_heads.size(), F_lens_and_origin_run, n, max_length_seen);
    if (domain != nullptr) {
        *domain = n;
    }
    if (max_length != nullptr) {
        *max_length = max_length_seen;
    }
    return {F_lens, F_images};
}
    
} // namespace orbit::rlbwt

#endif