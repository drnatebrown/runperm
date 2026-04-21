#ifndef _MOVE_SPLITTING_HH
#define _MOVE_SPLITTING_HH

#include "orbit/common.hpp"
#include "orbit/internal/ds/packed_vector.hpp"
#include <cmath>
#include <cassert>
#include <optional>
#include <iostream>

namespace orbit {

inline constexpr std::optional<double> DEFAULT_LENGTH_CAPPING = 8.0;
inline constexpr std::optional<ulint> DEFAULT_BALANCING = 16;

struct split_params {
    std::optional<double> length_capping;
    std::optional<ulint> balancing;

    split_params() : length_capping(DEFAULT_LENGTH_CAPPING), balancing(DEFAULT_BALANCING) {}
    split_params(std::optional<double> length_capping, std::optional<ulint> balancing)
    : length_capping(std::move(length_capping)), balancing(std::move(balancing)) {}

    bool operator==(const split_params& other) const {
        return length_capping == other.length_capping && balancing == other.balancing;
    }
    bool operator!=(const split_params& other) const { return !(*this == other); }

    size_t serialize(std::ostream& out) const {
        size_t written_bytes = 0;
        bool has_length_capping = length_capping.has_value();
        out.write(reinterpret_cast<const char*>(&has_length_capping), sizeof(has_length_capping));
        written_bytes += sizeof(has_length_capping);
        if (length_capping.has_value()) {
            double length_capping_value = length_capping.value();
            out.write(reinterpret_cast<const char*>(&length_capping_value), sizeof(length_capping_value));
            written_bytes += sizeof(length_capping_value);
        }
        bool has_balancing = balancing.has_value();
        out.write(reinterpret_cast<const char*>(&has_balancing), sizeof(has_balancing));
        written_bytes += sizeof(has_balancing);
        if (balancing.has_value()) {
            ulint balancing_value = balancing.value();
            out.write(reinterpret_cast<const char*>(&balancing_value), sizeof(balancing_value));
            written_bytes += sizeof(balancing_value);
        }
        return written_bytes;
    }

    void load(std::istream& in) {
        length_capping = std::nullopt;
        balancing = std::nullopt;
        bool has_length_capping;
        in.read(reinterpret_cast<char*>(&has_length_capping), sizeof(has_length_capping));
        if (has_length_capping) {
            double length_capping_value;
            in.read(reinterpret_cast<char*>(&length_capping_value), sizeof(length_capping_value));
            length_capping = length_capping_value;
        }
        bool has_balancing;
        in.read(reinterpret_cast<char*>(&has_balancing), sizeof(has_balancing));
        if (has_balancing) {
            ulint balancing_value;
            in.read(reinterpret_cast<char*>(&balancing_value), sizeof(balancing_value));
            balancing = balancing_value;
        }
    }
};

inline split_params NO_SPLITTING = split_params(std::nullopt, std::nullopt);
inline split_params DEFAULT_SPLITTING = split_params(DEFAULT_LENGTH_CAPPING, DEFAULT_BALANCING);
inline split_params ONLY_LENGTH_CAPPING = split_params(DEFAULT_LENGTH_CAPPING, std::nullopt);
inline split_params ONLY_BALANCING = split_params(std::nullopt, DEFAULT_BALANCING);

template<class int_vector_t>
struct split_result {
    int_vector_t lengths;
    int_vector_t img_rank_inv;
    ulint max_length;
};

class move_splitting {
public:

    template<class int_vector_t>
    inline static void split_by_length_capping(
        const int_vector_t& lengths, 
        const int_vector_t& img_rank_inv, 
        const ulint domain, 
        const double length_capping_factor, 
        split_result<int_vector_t>& result
    ) {
        assert(lengths.size() == img_rank_inv.size());
        assert(length_capping_factor > 0.0);

        double avg_run_length = static_cast<double>(domain) / static_cast<double>(lengths.size());
        ulint desired_max_allowed_length = static_cast<ulint>(std::ceil(avg_run_length * length_capping_factor));
        uchar bits = bit_width(desired_max_allowed_length);
        ulint max_allowed_length = max_val(bits);

        size_t new_intervals_upper_bound = std::ceil(static_cast<double>(lengths.size()) / length_capping_factor);

        split_by_max_allowed_length(lengths, img_rank_inv, domain, max_allowed_length, result, new_intervals_upper_bound);
    }

    template<class int_vector_t>
    inline static void split_by_max_allowed_length(
        const int_vector_t& lengths, 
        const int_vector_t& img_rank_inv,
        const ulint domain,
        const ulint max_allowed_length, 
        split_result<int_vector_t>& result,
        const size_t new_intervals_upper_bound = 0
    ) {
        assert(lengths.size() == img_rank_inv.size());
        assert(max_allowed_length > 0);
        size_t intervals_after_splitting_upper_bound = 
            (new_intervals_upper_bound == 0) ? domain - 1 : lengths.size() + new_intervals_upper_bound - 1;

        // For each original run, how many new intervals were added up to but not including this run?
        int_vector_t input_splits_exclusive_cumsum(lengths.size(), bit_width(new_intervals_upper_bound));
        
        // First pass to determine the number of intervals after splitting in input order, 
        // the number of new intervals added up to but not including each run, 
        // and the max length of the new intervals
        size_t num_intervals_after_splitting = 0;
        size_t cumulative_new_intervals = 0;
        result.max_length = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            input_splits_exclusive_cumsum[i] = cumulative_new_intervals;

            if (lengths[i] > max_allowed_length) {
                ulint remaining = lengths[i];
                while (remaining > 0) {
                    ulint chunk = std::min(remaining, max_allowed_length);
                    remaining -= chunk;
                    ++num_intervals_after_splitting;
                    ++cumulative_new_intervals;
                }
                --cumulative_new_intervals;
                result.max_length = max_allowed_length;
            } else {
                result.max_length = std::max(result.max_length, lengths[i]);
                ++num_intervals_after_splitting;
            }
        }

        result.lengths = int_vector_t(num_intervals_after_splitting, bit_width(result.max_length));
        result.img_rank_inv = int_vector_t(num_intervals_after_splitting, bit_width(num_intervals_after_splitting - 1));

        // Second pass to fill the lengths and img_rank_inv arrays, in output order
        size_t curr_img_rank_inv_idx = 0;
        for (size_t i = 0; i < lengths.size(); ++i) {
            size_t j = img_rank_inv[i];
            size_t length = lengths[j];
            size_t num_splits = 0;
            if (length > max_allowed_length) {
                ulint remaining = length;
                while (remaining > 0) {
                    ulint chunk = std::min(remaining, max_allowed_length);
                    remaining -= chunk;
                    result.lengths[j + input_splits_exclusive_cumsum[j] + num_splits] = chunk;
                    ++num_splits;
                }
                --num_splits;
            } else {
                result.lengths[j + input_splits_exclusive_cumsum[j]] = length;
            }
            
            // Fill the img_rank_inv array
            size_t curr_img_rank_inv_val = j + input_splits_exclusive_cumsum[j];
            for (size_t k = 0; k < num_splits + 1; ++k) {
                result.img_rank_inv[curr_img_rank_inv_idx] = curr_img_rank_inv_val + k;
                ++curr_img_rank_inv_idx;
            }
        }
    }

    inline static ulint balancing_upper_bound(ulint balancing_factor, ulint original_intervals) {
        return (((balancing_factor + 1) * original_intervals)/(balancing_factor - 1)) + 1;
    }

    template<class int_vector_t>
    inline static void split_by_balancing(const int_vector_t& lengths, const int_vector_t& img_rank_inv, const ulint domain, const ulint balancing_factor, split_result<int_vector_t>& result) {
        assert(lengths.size() == img_rank_inv.size());
        assert(balancing_factor >= 2);

        ulint intervals_upper_bound = balancing_upper_bound(balancing_factor, lengths.size());
        
        balance_state state;
        state.next_free_idx = lengths.size();
        state.balancing_factor = balancing_factor;
        state.balanced_up_to = 0;

        interval_list<int_vector_t> input_intervals(state, lengths.size(), intervals_upper_bound, domain);
        interval_list<int_vector_t> output_intervals(state, lengths.size(), intervals_upper_bound, domain);
        
        int_vector_t img_rank(lengths.size(), bit_width(lengths.size() - 1));
        initialize_lists(input_intervals, output_intervals, lengths, img_rank, img_rank_inv);
        
        // state.last_call_heavy = false;
        // state.skip_ahead_idx = 0;

        state.input_idx = 0; // index of input interval which intersects with the balanced_up_to_idx
        state.output_idx = 0; // index of output interval which intersects with the balanced_up_to_idx
        while (state.input_idx != input_intervals.END_IDX || state.output_idx != output_intervals.END_IDX) {
            ulint curr_input_start = input_intervals.get_start(state.input_idx);
            ulint curr_output_start = output_intervals.get_start(state.output_idx);
            ulint curr_input_length = (state.input_idx == input_intervals.END_IDX) ? 0 : input_intervals.get_start(input_intervals.get_next(state.input_idx)) - curr_input_start;
            ulint curr_output_length = (state.output_idx == output_intervals.END_IDX) ? 0 : output_intervals.get_start(output_intervals.get_next(state.output_idx)) - curr_output_start;

            // balance input
            if (curr_input_start < curr_output_start || state.output_idx == output_intervals.END_IDX) {
                std::optional<ulint> split_start = balance_input(input_intervals, output_intervals, state);
                update_balanced_up_to_idx(input_intervals, output_intervals, state, split_start, true);
            }
            // balance output
            else if (curr_output_start < curr_input_start || state.input_idx == input_intervals.END_IDX) {
                std::optional<ulint> split_start = balance_output(input_intervals, output_intervals, state);
                update_balanced_up_to_idx(input_intervals, output_intervals, state, split_start, false);
            }
            else {
                // balance input
                if (curr_input_length > curr_output_length) {
                    std::optional<ulint> split_start = balance_input(input_intervals, output_intervals, state);
                    update_balanced_up_to_idx(input_intervals, output_intervals, state, split_start, true);
                }
                // balance output
                else if (curr_output_length > curr_input_length) {
                    std::optional<ulint> split_start = balance_output(input_intervals, output_intervals, state);
                    update_balanced_up_to_idx(input_intervals, output_intervals, state, split_start, false);
                }
                else {
                    // No balancing needed
                    update_balanced_up_to_idx(input_intervals, output_intervals, state);
                }
            }
        }

        size_t num_intervals_after_splitting = state.next_free_idx;
        result.max_length = 0;

        ulint idx = 0;
        ulint curr_original_idx = 0;
        ulint cumulative_new_intervals = 0;
        int_vector_t input_splits_exclusive_cumsum(lengths.size(), bit_width(num_intervals_after_splitting - lengths.size()));
        for (size_t i = 0; i < num_intervals_after_splitting; ++i) {
            result.max_length = std::max(result.max_length, input_intervals.get_length(idx));
            if (idx < lengths.size()) {
                input_splits_exclusive_cumsum[curr_original_idx] = cumulative_new_intervals;
                ++curr_original_idx;
            }
            else {
                ++cumulative_new_intervals;
            }
            idx = input_intervals.get_next(idx);
        }
        assert(idx == input_intervals.END_IDX);

        result.lengths = int_vector_t(num_intervals_after_splitting, bit_width(result.max_length));
        result.img_rank_inv = int_vector_t(num_intervals_after_splitting, bit_width(num_intervals_after_splitting - 1));
        idx = 0;
        size_t steps_since_last_original = 0;
        ulint last_original_mapping = 0;
        for (size_t i = 0; i < num_intervals_after_splitting; ++i) {
            ulint curr_length = output_intervals.get_length(idx);
            if (idx < lengths.size()) {
                last_original_mapping = output_intervals.get_mapping(idx);
                steps_since_last_original = 0;
            }
            else {
                ++steps_since_last_original;
            }
            ulint curr_img_rank_inv = last_original_mapping + input_splits_exclusive_cumsum[last_original_mapping] + steps_since_last_original;
            result.img_rank_inv[i] = curr_img_rank_inv;
            result.lengths[curr_img_rank_inv] = curr_length;
            idx = output_intervals.get_next(idx);   
        }
        assert(idx == output_intervals.END_IDX);
    }

private:
    // TODO instead of END_IDX and domain, add to [END_IDX] with the start of domain to the input/output intervals
    struct balance_state {
        size_t next_free_idx;
        ulint balancing_factor;

        ulint balanced_up_to; // any value less than is balanced, less than or equal has pred values maintained
        size_t input_idx; // index of input interval which intersects with the balanced_up_to_idx
        size_t output_idx; // index of output interval which intersects with the balanced_up_to_idx

        // bool last_call_heavy; // whether the last call to balance was heavy (weight > 2 * balancing_factor)
        // ulint skip_ahead_idx; // where we scanned to when counting weight, don't need to here
    };

    DEFINE_ORBIT_COLUMNS(interval_cols, START, NEXT, PRED);

    template<class int_vector_t = int_vector>
    struct interval_list {
        packed_vector<interval_cols> orig_intervals;
        packed_vector<interval_cols> new_intervals;
        const int_vector_t* orig_mapping = nullptr;
        
        size_t domain;
        size_t END_IDX;

        interval_list(const balance_state& state, size_t original_intervals, size_t intervals_upper_bound, size_t domain) : domain(domain), END_IDX(intervals_upper_bound) {
            size_t new_intervals_upper_bound = intervals_upper_bound - original_intervals;
            uchar start_bits = bit_width(domain - 1); // start needs whole domain
            uchar orig_next_bits = bit_width(new_intervals_upper_bound); // next only points to new intervals, otherwise 0; need + 1 since 0 is reserved
            uchar pred_bits = bit_width(intervals_upper_bound - 1); // pred points over full range of intervals
            orig_intervals = packed_vector<interval_cols>(original_intervals, {start_bits, orig_next_bits, pred_bits});

            uchar new_next_bits = bit_width(intervals_upper_bound); // next points to next interval, or END_IDX if at end, over full range of intervals
            new_intervals = packed_vector<interval_cols>(new_intervals_upper_bound, {start_bits, new_next_bits, pred_bits});
            // we set orig_mapping later, since we can reuse the passed input img_rank_inv
        }

        inline void set_start(ulint idx, ulint start) {
            if (idx < orig_intervals.size()) {
                orig_intervals.set<interval_cols::START>(idx, start);
            }
            else {
                new_intervals.set<interval_cols::START>(idx - orig_intervals.size(), start);
            }
        }

        inline void init_next(ulint idx) {
            orig_intervals.set<interval_cols::NEXT>(idx, 0);
        }

        // Should only be used for updates, not initialization
        inline void set_next(ulint idx, ulint next_idx) {
            if (idx < orig_intervals.size()) {
                if (next_idx >= orig_intervals.size()) {
                    orig_intervals.set<interval_cols::NEXT>(idx, next_idx - orig_intervals.size() + 1);
                }
                else {
                    throw std::invalid_argument("Invalid next index");
                }
            }
            else {
                new_intervals.set<interval_cols::NEXT>(idx - orig_intervals.size(), next_idx);
            }
        }

        inline void set_pred(ulint idx, ulint pred_idx) {
            if (idx < orig_intervals.size()) {
                orig_intervals.set<interval_cols::PRED>(idx, pred_idx);
            }
            else {
                new_intervals.set<interval_cols::PRED>(idx - orig_intervals.size(), pred_idx);
            }
        }

        inline ulint get_start(ulint idx) {
            if (idx == END_IDX) {
                return domain;
            }
            else {
                if (idx < orig_intervals.size()) {
                    return orig_intervals.get<interval_cols::START>(idx);
                }
                else {
                    return new_intervals.get<interval_cols::START>(idx - orig_intervals.size());
                }
            }
        }

        inline ulint get_length(ulint idx) {
            if (idx == END_IDX) {
                return 0;
            }
            else {
                return get_start(get_next(idx)) - get_start(idx);
            }
        }

        inline ulint get_next(ulint idx) {
            if (idx < orig_intervals.size()) {
                ulint next_idx = orig_intervals.get<interval_cols::NEXT>(idx);
                // If next is 0, means next is just idx + 1
                if (next_idx == 0) {
                    if (idx == orig_intervals.size() - 1) {
                        return END_IDX;
                    }
                    else {
                        return idx + 1;
                    }
                }
                // Original is only non-zero if it points to a new interval, whose index is offset by the size of the original intervals
                else {
                    return orig_intervals.size() + next_idx - 1;
                }
            }
            // Otherwise, the new intervals hold the index of the next interval
            else {
                return new_intervals.get<interval_cols::NEXT>(idx - orig_intervals.size());
            }
        }

        inline ulint get_pred(ulint idx) {
            if (idx < orig_intervals.size()) {
                return orig_intervals.get<interval_cols::PRED>(idx);
            }
            else {
                return new_intervals.get<interval_cols::PRED>(idx - orig_intervals.size());
            }
        }

        inline ulint get_mapping(ulint idx) {
            if (idx < orig_intervals.size()) {
                return (*orig_mapping)[idx];
            }
            // If new interval, the mapping is always the index, since we added to input and output intervals in parallel
            else {
                return idx;
            }
        }
    };

    template<class int_vector_t>
    static inline void initialize_lists(interval_list<int_vector_t>& input_intervals, interval_list<int_vector_t>& output_intervals, const int_vector_t& lengths, int_vector_t& img_rank, const int_vector_t& img_rank_inv) {
        size_t curr_input_idx = 0;
        size_t curr_output_idx = 0;
        ulint curr_input_start = 0;
        ulint curr_output_start = 0;

        uchar mapping_bits = bit_width(lengths.size() - 1);
        output_intervals.orig_mapping = &img_rank_inv;
        input_intervals.orig_mapping = &img_rank;

        while (curr_input_idx < lengths.size() || curr_output_idx < lengths.size()) {
            ulint curr_input_length = (curr_input_idx == lengths.size()) ? 0 : lengths[curr_input_idx];
            ulint curr_output_length = (curr_output_idx == lengths.size()) ? 0 : lengths[img_rank_inv[curr_output_idx]];

            auto update_input_interval = [&](ulint output_pred) {
                input_intervals.set_start(curr_input_idx, curr_input_start);
                input_intervals.init_next(curr_input_idx);
                input_intervals.set_pred(curr_input_idx, output_pred);
            };
            auto update_output_interval = [&](ulint input_pred) {
                output_intervals.set_start(curr_output_idx, curr_output_start);
                output_intervals.init_next(curr_output_idx);
                output_intervals.set_pred(curr_output_idx, input_pred);
                
                img_rank[img_rank_inv[curr_output_idx]] = curr_output_idx;
            };
            
            if (curr_input_start < curr_output_start || curr_output_idx == lengths.size()) {
                update_input_interval(curr_output_idx - 1);
                curr_input_start += curr_input_length;
                curr_input_idx++;
            } else if (curr_output_start < curr_input_start || curr_input_idx == lengths.size()) {
                update_output_interval(curr_input_idx - 1);
                curr_output_start += curr_output_length;
                curr_output_idx++;
            } else {
                update_input_interval(curr_output_idx);
                update_output_interval(curr_input_idx);
                curr_input_start += curr_input_length;
                curr_output_start += curr_output_length;
                ++curr_input_idx;
                ++curr_output_idx;
            }
        }
    }

    template<class int_vector_t>
    static inline std::optional<ulint> balance_input(interval_list<int_vector_t>& input_intervals, interval_list<int_vector_t>& output_intervals, balance_state& state) {
        return balance(input_intervals, output_intervals, state.input_idx, state);
    }
    
    template<class int_vector_t>
    static inline std::optional<ulint> balance_output(interval_list<int_vector_t>& input_intervals, interval_list<int_vector_t>& output_intervals, balance_state& state) {
        return balance(output_intervals, input_intervals, state.output_idx, state);
    }
    
    // Does not set pred for the corresponding insertion into overlap_intervals
    template<class int_vector_t>
    static inline void split_interval(interval_list<int_vector_t>& to_balance_intervals, interval_list<int_vector_t>& overlap_intervals, ulint to_balance_idx, ulint overlap_idx, balance_state& state) {
        // Insert into to_balance_intervals
        ulint new_balance_interval_start = overlap_intervals.get_start(overlap_idx);
        to_balance_intervals.set_start(state.next_free_idx, new_balance_interval_start);
        to_balance_intervals.set_next(state.next_free_idx, to_balance_intervals.get_next(to_balance_idx));
        to_balance_intervals.set_pred(state.next_free_idx, overlap_idx);
    
        // Update existing intervals around the split interval
        to_balance_intervals.set_next(to_balance_idx, state.next_free_idx);
        overlap_intervals.set_pred(overlap_idx, state.next_free_idx);

        // Need to make corresponding insertion of mapped position to updated mapping values
        ulint delta = new_balance_interval_start - to_balance_intervals.get_start(to_balance_idx);
        ulint mapped_to_overlap_idx = to_balance_intervals.get_mapping(to_balance_idx);
        ulint new_overlap_idx = overlap_intervals.get_start(mapped_to_overlap_idx) + delta;

        // Insert into overlap_intervals
        overlap_intervals.set_start(state.next_free_idx, new_overlap_idx);
        overlap_intervals.set_next(state.next_free_idx, overlap_intervals.get_next(mapped_to_overlap_idx));

        // Update existing intervals around the mapped position
        overlap_intervals.set_next(mapped_to_overlap_idx, state.next_free_idx);

        ++state.next_free_idx;
    }
    
    // Returns the value of split in to_balance_intervals if a split was made, otherwise nullopt
    // TODO look at commented out skip ahead logic to avoid duplicate work, need to add first_call flag (also look at the state object)
    template<class int_vector_t>
    static inline std::optional<ulint> balance(interval_list<int_vector_t>& to_balance_intervals, interval_list<int_vector_t>& overlap_intervals, ulint to_balance_idx, balance_state& state) {
        ulint to_interval_start = to_balance_intervals.get_start(to_balance_idx);
        ulint to_interval_end = to_balance_intervals.get_start(to_balance_intervals.get_next(to_balance_idx));
        
        // ulint overlap_idx = (first_call && state.skip_ahead_idx != END_IDX) ? state.skip_ahead_idx : get_pred(to_balance_intervals, to_balance_idx);
        // ulint weight = (first_call && state.last_call_heavy) ? 0 : state.balancing_factor;
        ulint overlap_idx = to_balance_intervals.get_pred(to_balance_idx);
        ulint weight = 0;
        ulint split_candidate_idx = 0; // set to index of (balancing_factor + 1) position from start of overlap_idx
    
        // Find weight (number of starts in overlap_intervals that are contained in the interval to balance)
        ulint overlap_start = overlap_intervals.get_start(overlap_idx);
        while (overlap_start < to_interval_end && weight < 2 * state.balancing_factor) {
            if (overlap_start > to_interval_start) {
                ++weight;
            }
            
            // if (weight > 0 && first_call) {
            //     overlap_intervals.set<interval_cols::PRED>(overlap_idx, to_balance_idx);
            // }

            if (weight == state.balancing_factor + 1) {
                split_candidate_idx = overlap_idx;
            }

            overlap_idx = overlap_intervals.get_next(overlap_idx);
            overlap_start = overlap_intervals.get_start(overlap_idx);
        }
    
        // if (first_call) {
        //     state.skip_ahead_idx = overlap_idx;
        // }

        if (weight == 2 * state.balancing_factor) {
            split_interval(to_balance_intervals, overlap_intervals, to_balance_idx, split_candidate_idx, state);
            ulint new_overlap_idx = state.next_free_idx - 1;
            ulint new_to_balance_idx = state.next_free_idx - 1;
            ulint new_overlap_start = overlap_intervals.get_start(new_overlap_idx);
            ulint new_to_balance_start = to_balance_intervals.get_start(new_to_balance_idx);

            // Need to update pred values and make recursive call if the new overlap interval is in the balanced portion of the intervals
            if (new_overlap_start < state.balanced_up_to) {
                ulint mapped_to_overlap_idx = to_balance_intervals.get_mapping(to_balance_idx);
    
                // Gets the interval in to_balance_intervals that overlaps with the new overlap interval
                ulint new_overlap_pred_in_balance = overlap_intervals.get_pred(mapped_to_overlap_idx);
                ulint new_overlap_pred_in_balance_next = to_balance_intervals.get_next(new_overlap_pred_in_balance);
                while (to_balance_intervals.get_start(new_overlap_pred_in_balance_next) < new_overlap_start) {
                    new_overlap_pred_in_balance = new_overlap_pred_in_balance_next;
                    new_overlap_pred_in_balance_next = to_balance_intervals.get_next(new_overlap_pred_in_balance);
                }

                // Update predecessor value for new overlap interval
                overlap_intervals.set_pred(new_overlap_idx, new_overlap_pred_in_balance);

                // Update predecessor value for the interval in to_balance_intervals that overlaps with the new overlap interval
                ulint new_overlaps_in_balance = (to_balance_intervals.get_start(new_overlap_pred_in_balance) == new_overlap_start) ? new_overlap_pred_in_balance : new_overlap_pred_in_balance_next;
                // Get the end of the new overlap interval
                ulint new_overlap_end = overlap_intervals.get_start(overlap_intervals.get_next(new_overlap_idx));
                while (to_balance_intervals.get_start(new_overlaps_in_balance) < new_overlap_end) {
                    to_balance_intervals.set_pred(new_overlaps_in_balance, new_overlap_idx);
                    new_overlaps_in_balance = to_balance_intervals.get_next(new_overlaps_in_balance);
                }
                
                // balance(to_balance_intervals, overlap_intervals, new_overlap_pred_in_balance, state, false);
                balance(to_balance_intervals, overlap_intervals, new_overlap_pred_in_balance, state);
            }
            return new_to_balance_start;
        }
        return std::nullopt;
    }

    template<class int_vector_t>
    static inline void update_balanced_up_to_idx(interval_list<int_vector_t>& input_intervals, interval_list<int_vector_t>& output_intervals, balance_state& state, std::optional<ulint> split_start = std::nullopt, bool input_balance_step = true) {
        if (split_start.has_value()) {
            state.balanced_up_to = split_start.value();
        }
        else {
            if (input_balance_step) {
                state.balanced_up_to = input_intervals.get_start(input_intervals.get_next(state.input_idx));
            }
            else {
                state.balanced_up_to = output_intervals.get_start(output_intervals.get_next(state.output_idx));
            }
        }

        // if (input_balance_step) {
        //     state.input_idx = get_next(input_intervals, state.input_idx);
        // }
        // else {
        //     state.output_idx = get_next(output_intervals, state.output_idx);
        // }

        // update any predecessor values
        ulint prev_input_idx = state.input_idx;
        ulint prev_output_idx = state.output_idx;
        ulint curr_input_idx = (state.input_idx == input_intervals.END_IDX) ? input_intervals.END_IDX : input_intervals.get_next(state.input_idx);
        ulint curr_output_idx = (state.output_idx == output_intervals.END_IDX) ? output_intervals.END_IDX : output_intervals.get_next(state.output_idx);
        ulint curr_input_start = input_intervals.get_start(curr_input_idx);
        ulint curr_output_start = output_intervals.get_start(curr_output_idx);

        while ((curr_input_start <= state.balanced_up_to || curr_output_start <= state.balanced_up_to) && (curr_input_idx != input_intervals.END_IDX || curr_output_idx != output_intervals.END_IDX)) {
            if (curr_input_start < curr_output_start || curr_output_idx == output_intervals.END_IDX) {
                input_intervals.set_pred(curr_input_idx, prev_output_idx);
                state.input_idx = curr_input_idx;

                prev_input_idx = curr_input_idx;
                curr_input_idx = input_intervals.get_next(curr_input_idx);
                curr_input_start = input_intervals.get_start(curr_input_idx);   
            } else if (curr_output_start < curr_input_start || state.input_idx == input_intervals.END_IDX) {
                output_intervals.set_pred(curr_output_idx, prev_input_idx);
                state.output_idx = curr_output_idx;

                prev_output_idx = curr_output_idx;
                curr_output_idx = output_intervals.get_next(curr_output_idx);
                curr_output_start = output_intervals.get_start(curr_output_idx);
            } else {
                input_intervals.set_pred(curr_input_idx, curr_output_idx);
                output_intervals.set_pred(curr_output_idx, curr_input_idx);
                state.input_idx = curr_input_idx;
                state.output_idx = curr_output_idx;

                prev_input_idx = curr_input_idx;
                prev_output_idx = curr_output_idx;
                curr_input_idx = input_intervals.get_next(curr_input_idx);
                curr_output_idx = output_intervals.get_next(curr_output_idx);
                curr_input_start = input_intervals.get_start(curr_input_idx);
                curr_output_start = output_intervals.get_start(curr_output_idx);
            }
        }

        if (state.balanced_up_to == input_intervals.domain && curr_input_idx == input_intervals.END_IDX && curr_output_idx == output_intervals.END_IDX) {
            state.input_idx = input_intervals.END_IDX;
            state.output_idx = output_intervals.END_IDX;
        }
    }
};

} // namespace orbit
#endif