#!/usr/bin/env python

import scipy
import numpy
"""
window = [start, number_of_tags]
island = [current_chromosome_name, island_start, window_start + window_size, island_score, island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])
"""

# for input
# calculate_and_append_score_for_fdr(island_list_input, window_list_control_dict, input_lambda, window_size, NORMALIZATION_CONSTANT, 1)

# for control
#calculate_and_append_score_for_fdr(island_list_control, window_list_input_dict, control_lambda, window_size, 1, NORMALIZATION_CONSTANT)

calculate_and_append_score_for_fdr(island_list_input, window_list_control_dict, input_lambda, window_size, NORMALIZATION_CONSTANT, 1)
calculate_and_append_score_for_fdr(island_list_control, window_list_input_dict, control_lambda, window_size, 1, NORMALIZATION_CONSTANT)

"""
for island in island_list_input:
    while len(island) > 7:
        island.pop()
"""


def calculate_and_append_score_for_fdr(island_list_treatment, window_list_control_dict, lambdaa_treatment, window_size, normalization_4_treatment, normalization_4_control):
    """ Calculates score of island in treatment, compared with island in control at the same position"""

    previous_chr = ""
    for island in island_list_treatment:
        current_chr = island[0]
        island_start = island[1]
        island_end = island[2]
        island_reads_count = island[3]

        island_length = island_end - island_start
        windows_per_island = island_length / window_size

        island_lambda = island_reads_count * normalization_4_treatment
        lambdaa_treatment_normalized = lambdaa_treatment * windows_per_island * normalization_4_treatment
        # print('for 1')
        if current_chr != previous_chr:
            i = 0
            previous_chr = current_chr
        # count of control reads per treatment island
        control_reads_count = 0
        control_windows_per_island = 0

        while i < len(window_list_control_dict[current_chr]) - windows_per_island:
            window_list = window_list_control_dict[current_chr]
            # first window in island
            window_start = window_list[i][0]
            # print(last_window)

            if island_start > window_start:
                i += 1
                continue
            elif window_start > island_end:
                # if windows is over, check if all treatment_island was filled by windows
                # if not - cover the rest with average reads count (lambdaa_treatment)
                unfilled_windows = windows_per_island - control_windows_per_island
                control_reads_count += unfilled_windows * lambdaa_treatment_normalized  # lambda_control?
                break
            else:
                # number of control tags per island
                control_windows_per_island += 1
                number_of_reads_per_window = window_list[i][1]
                control_reads_count += number_of_reads_per_window
                i += 1

        control_lambda = control_reads_count * normalization_4_control
        lambdaa_treatment_counted = max(control_lambda, lambdaa_treatment_normalized)

        island.append(control_lambda)
        island.append(lambdaa_treatment_counted)
        island.append(lambdaa_treatment)


        # print(current_chr)
        # p = scipy.stats.poisson.cdf(island_lambda, lambdaa_treatment_counted)

        p = scipy.stats.poisson.ppf(island_lambda, lambdaa_treatment_counted)
        island.append(p)
        # sometimes p == nan
        if p <= 1e-320:
            score_for_fdr = 1000
        else:
            score_for_fdr = -numpy.log(p)

        island.append(score_for_fdr)

        # compare score_for_fdr with log(p0):
        #if score_for_fdr > -numpy.log(0.1):


# write_output(outfile, island_list_input, ISLAND_SCORE_THRESHOLD)

def make_scores_dict(island_list):
    # {score : [island_position_in_list]}
    scores_dict = {}
    for i in range(len(island_list)):
        island = island_list[i]
        score_for_fdr = island[10]

        island_position_in_list = i
        # chr_name = island[0]

        try:
            scores_dict[score_for_fdr].append(island_position_in_list)
        except KeyError:
            scores_dict[score_for_fdr] = []
            scores_dict[score_for_fdr].append(island_position_in_list)
    return scores_dict

# calculate_and_append_fdr(island_list_input, island_list_control)


def calculate_and_append_fdr(island_list_treatment, island_list_control):
    treatment_scores_dict = make_scores_dict(island_list_treatment)
    treatment_scores = treatment_scores_dict.keys()
    treatment_scores.sort(reverse=True)

    for island in island_list_treatment:
        island.append("NA")

    control_scores = []
    for island in island_list_control:
        score_for_fdr = island[10]
        control_scores.append(score_for_fdr)
    control_scores.sort(reverse=True)

    negative_islands = 0
    true_islands = 0

    for t_score in treatment_scores:
        true_islands += len(treatment_scores_dict[t_score])
        while negative_islands < len(control_scores) and t_score <= control_scores[negative_islands]:
            negative_islands += 1
        fdr = 100.0 * negative_islands / true_islands

        for position in treatment_scores_dict[t_score]:
            if fdr > 100:
                fdr = 100
            island_list_treatment[position][11] = fdr


# for control
# write_output('/media/user/DISK1/SICER_project/BAM_files/our_control/with_control_and_FDR_calculate_fdr_CONTROL_score_ppf.bed', island_list_control, ISLAND_SCORE_THRESHOLD)

"""
for i in island_list_input:
    while len(i) > 7:
        i.pop()


control_scores = set()

for i in island_list_control:
    control_scores.add(i[7])

input_scores = set()
for i in island_list_input:
    input_scores.add(i[7])
"""