#!/usr/bin/env python

import logging
import scipy
import scipy.stats

import pre_counting
import islands


def broadpeaks_wo_control(bam_path, WINDOW_SIZE, GAP, EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0):

    chromosomes_info = pre_counting.get_chromosomes_info(bam_path)

    logging.info("\nStep 1 of 4\nCOUNTING UNIQUE READS\n")
    input_unique_reads_count = pre_counting.count_unique_reads(bam_path, chromosomes_info)

    # Effective genome length (L)
    effective_length = pre_counting.count_effective_length(EFFECTIVE_PROPORTION, chromosomes_info)

    # Lambda for poisson distribution
    lambdaa = pre_counting.count_lambda(input_unique_reads_count, WINDOW_SIZE, effective_length)

    # Minimum #reads in a window for eligibility
    # Formula (1), finding l0
    l0 = scipy.stats.poisson.ppf(1 - p0, lambdaa)
    NORMALIZATION_CONSTANT = 1
    logging.info("\nWindow read threshold is {} reads, \ni.e. {} is minimum number of reads in window "
                 "to consider this window `eligible` with Poisson distribution p-value {}".format(l0, l0, p0))

    logging.info("\nStep 2 of 4\nMAKING WINDOW LIST\n")
    (window_list, window_list_dict) = islands.make_windows_list(bam_path, chromosomes_info, l0, WINDOW_SIZE, GAP,
                                            input_unique_reads_count, NORMALIZATION_CONSTANT)

    logging.info("\nStep 3 of 4\nMAKING ISLAND LIST\n")
    island_list = islands.make_islands_list(window_list, lambdaa, WINDOW_SIZE, l0, chromosomes_info,
                                            ISLAND_SCORE_THRESHOLD)

    return island_list
