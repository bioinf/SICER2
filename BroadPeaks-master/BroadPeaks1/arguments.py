#!/usr/bin/env python

import sys
import logging
import os.path

import colorized_log


def check_input(input_path):
    """
    :param input_path: (str)
    :return: input_path
    """
    # Check if there is a file at entered dir
    if not os.path.isfile(input_path):
        sys.exit("No BAM file specified in input '{}' or there is no such a file".format(input_path))
    # Check if it is a BAM file
    if input_path[-4:] != '.bam':
        sys.exit("`{}` is other file type, not BAM. This tool works only with BAM-files as input (*.bam)".
                 format(input_path))
    return input_path


def make_log(input_path, log_dir, log_name, merge_log, stop_verbosity, window_size, gap_size):
    # log_path = input_path[:-4]
    # specifying directory and name of LOG file
    extension = ".log"
    if log_name:
        if os.path.isdir(log_dir):
            log_path = log_dir + "/" + log_name + extension
        else:
            log_path = os.path.dirname(input_path) + "/" + log_name + extension
    else:
        prefix = "_W{}_G{}_log".format(window_size, gap_size * window_size)
        if os.path.isdir(log_dir):
            log_path = log_dir + "/" + os.path.basename(input_path)[:-4] + prefix + extension
        else:
            log_path = os.path.dirname(input_path) + "/" + os.path.basename(input_path)[:-4] + prefix + extension
    # making log file
    if merge_log:
        logging.basicConfig(filename=log_path, level=logging.DEBUG,
                            format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
    else:
        logging.basicConfig(filename=log_path, level=logging.DEBUG,
                            format='%(asctime)s : %(levelname)s : %(message)s', datefmt='%m/%d/%Y %H:%M:%S',
                            filemode='w')
    # deciding whether to print logs to console
    if not stop_verbosity:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
    # formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # console.setFormatter(formatter)
        formatter = colorized_log.ColoredFormatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
    print ("\nLOG file for this run is here: \n{}\n".format(log_path))


def check_control(control_path):
    if control_path != "":
        if not os.path.isfile(control_path):
            logging.error("No control BAM file specified in input '{}' or there is no such a file".format(control_path))
            sys.exit("`{}` is not a path to BAM file. \n More information in LOG file for this run".format(control_path))
        if control_path[-4:] != '.bam':
            logging.error("`{}` is other file type, not BAM. This tool works only with BAM-files as control input (*.bam)".
                          format(control_path))
            sys.exit("`{}` is not a BAM file. \n More information in `{}`".format(control_path))
    return control_path


def check_p_value(p0):
    if p0 <= 0 or p0 >= 1:
        logging.error("`{}` is incorrect p-value. p-value has to be in range(0.0, 1.0)".format(p0))
        sys.exit()
    return p0


def check_effective_proportion(effective_proportion):
    if effective_proportion <= 0 or effective_proportion > 1:
        logging.error("`{}` is incorrect proportion of effective genome length. \n "
                      "Proportion of effective genome length has to be in range(0.0, 1.0)".format(effective_proportion))
        sys.exit()
    return effective_proportion


def check_outfile(output_dir, output_name, input_path, window_size, gap_size):
    if not output_dir:
        output_dir = os.path.dirname(input_path)
    if not output_name:
        output_name = os.path.basename(input_path)[:-4] + "_W{}_G{}".format(window_size, window_size * gap_size)
    # must test validity of output_name as filename
    outfile = output_dir + "/" + output_name + "_peaks.bed"
    return outfile


def write_run_information(input_path, window_size, gap_size, island_score_threshold, effective_proportion,
                          control_path, outfile):
    logging_string = "You have started tool for ChiP-seq analysis to find broad peaks with following parameters:" \
                     "\n\n1. Input file : {}".format(input_path)
    logging_string += "\n2. Window size :" + str(window_size) + " bp"
    if window_size == 200:
        logging_string += " (default)"
    logging_string += "\n3. Gap size : {} bp".format(window_size * gap_size)
    if gap_size == 1:
        logging_string += " (default)"
    logging_string += "\n4. Island score threshold : {}".format(island_score_threshold)
    if island_score_threshold == 0:
        logging_string += " (default)"
    logging_string += "\n5. Proportion of effective genome length : {}".format(effective_proportion)
    if effective_proportion == 0.77:
        logging_string += " (default)"
    if control_path:
        logging_string += "\n6. Control file : {}".format(control_path)
    else:
        logging_string += "\n6. Control file : not specified (default)"
    logging_string += "\n7. Output file : {}".format(outfile)
    logging.info(logging_string + "\n")
