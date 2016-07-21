#!/usr/bin/env python

import logging
import scipy
import scipy.stats

import pre_counting
import islands
import FDR_calculation


def broadpeaks_with_control(bam_path, control_path, window_size, gap, EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0):

    chromosomes_info = pre_counting.get_chromosomes_info(bam_path)
    control_chromosomes_info = pre_counting.get_chromosomes_info(control_path)

    logging.info("\nStep 1 of 4\nCOUNTING UNIQUE READS\n")
    logging.info("\nFor input file\n")
    # input_unique_reads_count = 6300518
    input_unique_reads_count = pre_counting.count_unique_reads(bam_path, chromosomes_info)
    logging.info("\nFor control file\n")
    # control_unique_reads_count = 3758349
    control_unique_reads_count = pre_counting.count_unique_reads(control_path, control_chromosomes_info)

    # Effective genome length (L)
    # effective_length = 2383684366.91
    effective_length = pre_counting.count_effective_length(EFFECTIVE_PROPORTION, chromosomes_info)

    # Lambda for poisson distribution
    # input_lambda = 0.53
    input_lambda = pre_counting.count_lambda(input_unique_reads_count, window_size, effective_length)
    # control_lambda = 0.32
    control_lambda = pre_counting.count_lambda(control_unique_reads_count, window_size, effective_length)

    # Minimum #reads in a window for eligibility
    # Formula (1), finding l0
    input_l0 = scipy.stats.poisson.ppf(1 - p0, input_lambda)
    control_l0 = scipy.stats.poisson.ppf(1 - p0, control_lambda)
    logging.info("\nWindow read threshold is {} reads, \ni.e. {} is minimum number of reads in window "
                 "to consider this window `eligible` with Poisson distribution p-value {}".format(input_l0, input_l0, p0))


    logging.info("\nStep 2 of 4\nMAKING WINDOW LIST\n")

    # for two libraries are independent, we do not scale them here
    NORMALIZATION_CONSTANT = 1
    logging.info("\nFor input file\n")
    (window_list_input, window_list_input_dict) = islands.make_windows_list(bam_path, chromosomes_info, input_l0, window_size, gap, input_unique_reads_count, NORMALIZATION_CONSTANT)
    logging.info("\nFor control file\n")
    (window_list_control, window_list_control_dict) = islands.make_windows_list(control_path, chromosomes_info, control_l0, window_size, gap, control_unique_reads_count, NORMALIZATION_CONSTANT)

    #window_list = islands.modify_window_list_based_on_control(control_path, chromosomes_info, l0, window_size, gap, input_unique_reads_count, control_unique_reads_count, window_list_temp)


    logging.info("\nStep 3 of 4\nMAKING ISLAND LIST\n")

    logging.info("\nFor input file\n")
    island_list_input = islands.make_islands_list(window_list_input, input_lambda, window_size, input_l0, chromosomes_info, ISLAND_SCORE_THRESHOLD)
    logging.info("\nFor control file\n")
    island_list_control = islands.make_islands_list(window_list_control, control_lambda, window_size, control_l0, control_chromosomes_info, ISLAND_SCORE_THRESHOLD)

    """
    # with switching tracks
    calculate_fdr(island_list_input, window_list_control)
    calculate_fdr(island_list_control, window_list_input)


    island_list = islands.find_unintersected_islands(island_list_input,island_list_control)

    # calculate FDR
    FDR = (len(island_list_control) - (len(island_list_input)-len(island_list)))/len(island_list)
    logging.info("\nFDR is {} reads, \n".format(FDR))
    """
    # appending FDR to island_list_input
    # normalization to smaller dataset
    if input_unique_reads_count >= control_unique_reads_count:
        NORMALIZATION_CONSTANT = float(control_unique_reads_count) / input_unique_reads_count
        FDR_calculation.calculate_and_append_score_for_fdr(island_list_input, window_list_control_dict, input_lambda, window_size, NORMALIZATION_CONSTANT, 1)
        FDR_calculation.calculate_and_append_score_for_fdr(island_list_control, window_list_input_dict, control_lambda, window_size, 1, NORMALIZATION_CONSTANT)
    else:
        NORMALIZATION_CONSTANT = float(input_unique_reads_count) / control_unique_reads_count
        FDR_calculation.calculate_and_append_score_for_fdr(island_list_input, window_list_control_dict, input_lambda, window_size, 1, NORMALIZATION_CONSTANT)
        FDR_calculation.calculate_and_append_score_for_fdr(island_list_control, window_list_input_dict, control_lambda, window_size, NORMALIZATION_CONSTANT, 1)

    FDR_calculation.calculate_and_append_fdr(island_list_input, island_list_control, 7)
    FDR_calculation.calculate_and_append_fdr(island_list_input, island_list_control, 8)

    return island_list_input
