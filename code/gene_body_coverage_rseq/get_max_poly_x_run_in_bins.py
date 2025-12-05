#here we calculate longest poly A run for every gene
import os
import numpy as np
import argparse, sys

def calculate_poly_x_regions(sequence, nucleotides = 'A', allow_mismatches_without_breaking_streak = 1):
    """
    Calculate the regions of consecutive nucleotides in a given sequence.

    Parameters:
    sequence (str): The input nucleotide sequence.
    nucleotides (str): The nucleotide(s) to look for. Default is 'A'.
    allow_mismatches_without_breaking_streak (int): Number of allowed mismatches within a streak. Default is 1.

    Returns:
    list of tuples: Each tuple contains the start and end indices of a region.
    """
    

    regions = []
    streak_lengths = []
    start = None
    mismatch_count = 0
    streak_length = 0

    for i, nucleotide in enumerate(sequence):
        if nucleotide in nucleotides:
            if start is None:
                start = i
            mismatch_count = 0
            streak_length += 1
        else:
            if start is not None:
                mismatch_count += 1
                if mismatch_count > allow_mismatches_without_breaking_streak:
                    regions.append((start, i - mismatch_count))
                    streak_lengths.append(streak_length)
                    start = None
                    mismatch_count = 0
                    streak_length = 0

    if start is not None:
        regions.append((start, len(sequence) - 1))
        streak_lengths.append(streak_length)
        
    return regions, streak_lengths

def get_bin_ends(l, n_bins = 100):
    """
    Return a list of positions (1-based integers) dividing a sequence of length l
    into the same 100 percentile bins used by the original code.
    If l < 100, returns list(range(1, l+1)).
    """
    if l <= 0:
        return []
    if l < n_bins:
        return list(range(1, l+1))
    bins = []
    for i in range(1, n_bins+1):
        k = l * (i / n_bins)
        pos = int(round(k)) - 1  # matches original interpolation + int(round(...))
        bins.append(pos)
    return bins

def assign_longest_streak_to_bin(bins, regions, streak_lengths):
    
    region_lower_bounds = np.array([r[0] for r in regions])
    region_upper_bounds = np.array([r[1] for r in regions])
    max_streaks_in_bins = []
    #print(bins)
    #print(region_lower_bounds)
    #print(region_upper_bounds)
    for i, _ in enumerate(bins):
        if i > 0:
            bin_lower_pos = bins[i-1]
        else:
            bin_lower_pos = 0
        bin_upper_pos = bins[i]
        
        if np.max((region_upper_bounds >= bin_lower_pos) & (region_lower_bounds <= bin_upper_pos)) == 1:
            lowest_region_index = np.argmax((region_upper_bounds >= bin_lower_pos) & (region_lower_bounds <= bin_upper_pos))
        else:
            lowest_region_index = 'nan' #for debugging
            max_streaks_in_bins.append(0)
            continue
        furthest_region_index = np.sum((region_lower_bounds <= bin_upper_pos)) - 1
        #print(f'i : {i}')
        #print(f'bin_lower_pos : {bin_lower_pos}')
        #print(f'bin_upper_pos : {bin_upper_pos}')
        #print(f'lowest_region_index : {lowest_region_index}')
        #print(f'furthest_region_index : {furthest_region_index}')
        
        max_streaks_in_bins.append(int(np.max(streak_lengths[lowest_region_index:(furthest_region_index + 1)])))

    return max_streaks_in_bins


def get_max_poly_x_run_in_bins(sequence, nucleotides = 'A', allow_mismatches_without_breaking_streak = 1, n_bins = 100):
        regions, streak_lengths = calculate_poly_x_regions(sequence, nucleotides = nucleotides, allow_mismatches_without_breaking_streak = allow_mismatches_without_breaking_streak)
        bins = get_bin_ends(len(sequence), n_bins = n_bins)
        return assign_longest_streak_to_bin(bins, regions, streak_lengths)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate longest poly X run in each bin for input sequence.')
    parser.add_argument('-s', '--sequence', type=str, help='The sequence we analyze.')
    parser.add_argument('-n', '--nucleotides', type=str, default='A', help='The nucleotides which the streaks we look for are composed of. Default is "A".')
    parser.add_argument('-m', '--mismatches_allowed', type=int, default=1, help='Number of allowed mismatches within a streak. Default is 1.')
    parser.add_argument('-b', '--bin_number', type=int, default=100, help='Number of bins to divide the sequence into. Default is 100.')
    
    args = parser.parse_args()
    if args.sequence:
        sequence = args.sequence.strip()
    else:
        sequence = sys.stdin.read().strip()
    if not sequence:
        parser.error("No sequence provided via -s/--sequence or stdin.")
        
    nucleotides = args.nucleotides
    mismatches_allowed = args.mismatches_allowed
    bin_number = args.bin_number
    max_streaks = get_max_poly_x_run_in_bins(sequence, nucleotides=nucleotides, allow_mismatches_without_breaking_streak=mismatches_allowed, n_bins=bin_number)
    
    print(f"Longest streaks made up of '{nucleotides}' nucleotides per bin for {bin_number} bins ({mismatches_allowed} mismatches allowed without breaking streak): ")
    print(max_streaks)