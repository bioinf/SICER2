#!/usr/bin/env python

import argparse
import time
import scipy.stats
import numpy
import pre_counting

import arguments
import islands
import logging

startTime = time.time()

# ALL TO main()

parser = argparse.ArgumentParser(description="Tool for ChiP-seq analysis to find broad peaks")

parser.add_argument('infile', help="Path to `input.bam` file", type=str)
parser.add_argument('-w', dest='window_size', help="Window size (bp).  DEFAULT: 200", type=int, default=200)
parser.add_argument('-g', dest='gap', help="Gap size shows how many windows could be skipped. DEFAULT: 1",
                    type=int, default=1, choices=[1, 2, 3])
parser.add_argument('-t', dest='threshold', help="Island score threshold. DEFAULT: 0", type=int, default=0)
parser.add_argument('-o', dest='outdir', help="Path to directory for output `*_peaks.bed` file. "
                                              "DEFAULT: output will be in the same directory as `input.bam`",
                    type=str, default="")
parser.add_argument('-n', dest="output_name", help="Specify output name. "
                                                   "DEFAULT : an input file name + "
                                                   "`W<window size (bp)>_G<gap size (bp)>_peaks.bed`", type=str)
parser.add_argument('-e', help="Proportion of effective genome length; has to be in range(0.0, 1.0) DEFAULT: 0.77",
                    type=float, default=0.77)
parser.add_argument('-c', dest='control', help="Path to `control.bam` file. DEFAULT: no control file",
                    type=str, default="")
parser.add_argument('-log_name', help="Specify LOG file name."
                                       "DEFAULT : `input_log.log`", type=str, default="")
parser.add_argument('-log_dir', help="Specify path to directory where to write LOG file."
                                       "DEFAULT : log file will be in the same directory as output",
                    type=str, default="")
parser.add_argument('--merge_log', help="Merge logs from all runs in one LOG file. "
                                  "DEFAULT : LOG file contains information only from the last run", action='store_true')
parser.add_argument('--stop_verbosity', help="Stops printing logs to terminal, just to LOG file", action='store_true')


args = parser.parse_args(['/media/user/DISK1/SICER_project/BAM_files/H3K4Me3_test.bam', '-n', 'with_control_and_FDR_calculate_fdr', '-c', '/media/user/DISK1/SICER_project/BAM_files/test_control.bam', '-o', '/media/user/DISK1/SICER_project/BAM_files/our_control'])

bam_path = arguments.check_input(args.infile)
WINDOW_SIZE = args.window_size
GAP = args.gap
if args.outdir:
    arguments.make_log(bam_path, args.outdir, args.log_name, args.merge_log, args.stop_verbosity, WINDOW_SIZE, GAP)
else:
    arguments.make_log(bam_path, args.log_dir, args.log_name, args.merge_log, args.stop_verbosity, WINDOW_SIZE, GAP)
EFFECTIVE_PROPORTION = arguments.check_effective_proportion(args.e)
ISLAND_SCORE_THRESHOLD = args.threshold
outfile = arguments.check_outfile(args.outdir, args.output_name, bam_path, WINDOW_SIZE, GAP)
control_path = arguments.check_control(args.control)
# p0 = arguments.check_p_value(args.p_value)
p0 = 0.1


chromosomes_info = pre_counting.get_chromosomes_info(bam_path)
control_chromosomes_info = pre_counting.get_chromosomes_info(control_path)

window_size = 200
gap = 1
p0 = 0.1

effective_length = 2383684366.91
control_unique_reads_count = 3758349

input_lambda = 0.53
control_lambda = 0.32
input_l0 = scipy.stats.poisson.ppf(1 - p0, input_lambda)
control_l0 = scipy.stats.poisson.ppf(1 - p0, control_lambda)

input_unique_reads_count = 6300518


NORMALIZATION_CONSTANT = 1
(window_list_input, window_list_input_dict) = islands.make_windows_list(bam_path, chromosomes_info, input_l0, window_size, gap, input_unique_reads_count, NORMALIZATION_CONSTANT)
(window_list_control, window_list_control_dict) = islands.make_windows_list(control_path, chromosomes_info, control_l0, window_size, gap, control_unique_reads_count, NORMALIZATION_CONSTANT)


island_list_input = islands.make_islands_list(window_list_input, input_lambda, window_size, input_l0, chromosomes_info, ISLAND_SCORE_THRESHOLD)
island_list_control = islands.make_islands_list(window_list_control, control_lambda, window_size, control_l0, control_chromosomes_info, ISLAND_SCORE_THRESHOLD)

NORMALIZATION_CONSTANT = float(control_unique_reads_count) / input_unique_reads_count
