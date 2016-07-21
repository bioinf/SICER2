#!/usr/bin/env python

import numpy
import scipy
import pysam
import logging

# Makes a simple list of windows, where each window is a list [WindowStart, ReadsInWindow].
# Chromosomes are separated by [-1,-1] window
# Sequences of ineligible windows longer than GAP+1 are not stored


def make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, unique_reads_count, normalization_coef):
    logging.info("Making eligible windows of {} bp with allowed gap_size {} bp".format(window_size, window_size*gap))
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    window_list = []
    window_list_dict = {}
    # logging.info("chromosome_name, chromosome_size, total_number_of_eligible_windows_on_chromosome")
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        previous_read_strand = 0
        current_chromosome_name = chromosome[0]
        window_list_dict[current_chromosome_name] = []
        add_window = window_list_dict[current_chromosome_name].append
        pop_window = window_list_dict[current_chromosome_name].pop
        current_chromosome_size = int(chromosome[1])
        # print([current_chromosome_name, current_chromosome_size, len(window_list)])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        gap_count = 0
        window_start = 0
        window_reads_count = 0
        # just for precise number of windows counting
        chr_window = 0
        i = 0
        for read in all_reads_in_chromosome:
            # print((list(all_reads_in_chromosome)))
            read_str = str(read)
            # read strand: 0 = +         16 = -
            read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            # filtering redundant reads
            if (beginning_of_the_read != beginning_of_the_previous_read) or (read_strand != previous_read_strand):
                beginning_of_the_previous_read = beginning_of_the_read
                previous_read_strand = read_strand
                # Ground state: gap_count <= GAP
                gap_flag = True
                while True:
                    # if read in window
                    #print(window_start)
                    if window_start <= beginning_of_the_read < window_start + window_size:
                        window_reads_count += 1
                         # print(window_reads_count)
                        break
                    # if read before window: NEVER ENTERING THIS CONDITION
                    elif beginning_of_the_read < window_start:
                        break
                    else:
                        window_reads_count = int(float(window_reads_count)*normalization_coef)
                        if window_reads_count < l0:
                            gap_count += 1
                        else:
                            gap_flag = False
                            gap_count = 0
                        # * unique_reads_count/1000000 is for normalization per million reads
                        # now we are able to compare control and sample

                        # NEED TO CHANGE LAMBDA AND N
                        window_list.append([window_start, window_reads_count])
                        add_window([window_start, window_reads_count])
                        chr_window += 1
                        # / (unique_reads_count/1000000)])
                        # If we have a g+1 sized GAP, go and delete last g windows
                        #print(window_list)
                        # print(window_list_dict)
                        if gap_count > gap or gap_flag:
                            gap_flag = True
                            while gap_count > 0:
                                window_list.pop()
                                pop_window()
                                gap_count -= 1
                                chr_window -= 1
                        window_start += window_size
                        window_reads_count = 0
        # appends island if this is the last window
        if window_reads_count != 0:
            window_list.append([window_start, window_reads_count])
            add_window([window_start, window_reads_count])
        # Next chromosome marker just in case
        window_list.append([-1, -1])
        i += 1
        logging.info("On {} there are {} eligible windows".format(current_chromosome_name, chr_window))
    # candidate windows == eligible + ineligible(below gap)
    logging.info("\nThere are {} candidate windows".format(len(window_list) - i))
    window_list.append([1, 1])
    bamfile.close()
    # print(window_list)
    # print(window_list_dict)
    return (window_list, window_list_dict)


"""
# unused function for the other method of accounting for control
def modify_window_list_based_on_control(control_path, chromosomes_info, l0, window_size, gap, unique_reads_count, control_unique_reads_count, window_list_wo_control):
    logging.info("Making window list based on control")
    bamfile = pysam.AlignmentFile(control_path, 'rb')
    NORMALIZATION_CONSTANT = float(unique_reads_count)/float(control_unique_reads_count)
    window_list = []
    i = 0
    # logging.info("chromosome_name, chromosome_size, total_number_of_eligible_windows_on_chromosome")
    for chromosome in chromosomes_info:
        beginning_of_the_previous_read = 0
        previous_read_strand = 0
        current_chromosome_name = chromosome[0]
        current_chromosome_size = int(chromosome[1])
        all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
        window_start = window_list_wo_control[i][0]
        control_window_reads_count = 0

        for read in all_reads_in_chromosome:
            read_str = str(read)
            # read strand: 0 = +         16 = -
            read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
            beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
            # filtering redundant reads
            if (beginning_of_the_read != beginning_of_the_previous_read) or (read_strand != previous_read_strand):
                beginning_of_the_previous_read = beginning_of_the_read
                previous_read_strand = read_strand

                if window_start <= beginning_of_the_read and beginning_of_the_read < window_start + window_size:
                    control_window_reads_count += 1
                    break
                elif beginning_of_the_read < window_start:
                    break
                else:
                    window_reads_count = window_list_wo_control[i][1] - int(control_window_reads_count*NORMALIZATION_CONSTANT)
                    window_list.append([window_start, window_reads_count])
                    i += 1
                    window_start = window_list_wo_control[i][0]
                    #next chromosome check
                    if window_start == -1:
                        window_list.append([-1,-1])
                        window_start = window_list_wo_control[i][0]
                        i += 1
                        break
        #go and do next chrom

    #prune new window list of gaps
    gap_count = 0
    gap_flag = True
    window_list_new = []
    for window in window_list:
        window_list_new.append(window)
        if window[1] < l0:
            gap_count += 1
        else:
            gap_flag = False
        if gap_count > gap or gap_flag:
            gap_flag = True
            while gap_count > 0:
                window_list_new.pop()
                gap_count -= 1

    return window_list_new
"""



def calculate_window_score(reads_in_window, lambdaa, l0):
    # sometimes 0 and therefore inf in -log  is generated
    if reads_in_window >= l0:
        temp = scipy.stats.poisson.pmf(reads_in_window, lambdaa)
        if temp < 1e-320:
            window_score = 1000
        else:
            window_score = -numpy.log(temp)
    else:
        window_score = 0
    return window_score


def make_islands_list(window_list, lambdaa, window_size, l0, chromosomes_info, island_score_threshold):
    chromosome_counter = 0
    current_chromosome_name = chromosomes_info[chromosome_counter][0]
    islands_list = []
    island_score = 0
    island_number_of_reads = 0
    island_number_of_gaps = 0
    window_start = window_list[0][0] - window_size
    island_start = window_list[0][0]
    #zero_score_islands = 0
    for i, window in enumerate(window_list):
        # i == # in list, window == [window_start_position, number_of_reads_per_window]
        window_start_new = window[0]
        # what is window[1]?
        number_of_reads = window[1]

        # New chromosome check: [-1  -1] window separates chomosomes.
        if window_start_new == -1:
            # append the previous island
            island_number_of_reads += window_list[i-1][1]
            island_score += calculate_window_score(window_list[i-1][1], lambdaa, l0)
            # sometimes -200 appeared
            if island_start < 0:
                island_start = 0
            islands_list.append([current_chromosome_name, island_start, window_start + window_size, island_score,
                                 island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])

            # reset parameters for the new island
            chromosome_counter += 1
            # switch the chromosome name to next one
            if chromosome_counter < len(chromosomes_info):
                current_chromosome_name = chromosomes_info[chromosome_counter][0]
                window_start = window_list[i + 1][0] - window_size
                island_start = window_list[i + 1][0] - window_size
                island_score = 0
                island_number_of_reads = 0
                island_number_of_gaps = 0
            # if all chromosomes were observed
            elif chromosome_counter == len(chromosomes_info):
                break

        else:
            # if the next window belongs to the next island:
            if window_start_new != window_start + window_size:
                # Special case for the one-window island:
                if window_start == island_start:
                    island_score = calculate_window_score(number_of_reads, lambdaa, l0)
                    island_number_of_reads = number_of_reads
                island_length = (window_start + window_size - island_start)/window_size
                if island_start < 0:
                    island_start = 0
                islands_list.append([current_chromosome_name, island_start,
                                    window_start + window_size, island_score, island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])
                island_score = 0
                island_number_of_reads = 0
                island_number_of_gaps = 0
                island_start = window_start_new
            else:
                # Set score
                window_score = calculate_window_score(number_of_reads, lambdaa, l0)
                if number_of_reads < l0:
                    island_number_of_gaps += 1
                island_score += window_score
                island_number_of_reads += number_of_reads
            window_start = window_start_new
    # islands_list.append([current_chromosome_name, island_start, window_start + window_size, island_score,
                         #island_number_of_reads, island_length - island_number_of_gaps, island_number_of_gaps])
    logging.info("There are {} islands found".format(len(islands_list)))
    #print(zero_score_islands)
    # print(islands_list)
    return islands_list


"""
# diff function returns difference between two island lists, ie islands that are in the first list but not in the second.
# both lists are assumed to be sorted by chromosome and by island_start
def diff(island_list, island_list_2):
    i = 0
    final_islands = []
    # iterate through 2nd list, which is a subset of first
    for island in island_list_2:
        # if islands don't match, final result will contain the island from the first list
        while i< len(island_list):
            if island_list[i] != island:
                final_islands.append(island_list[i])
                i+=1
            else:
                i+=1
                break
    # some islands from the first list may be left after the iteration through second list is over. we append them to the result as well.
    while i< len(island_list):
            final_islands.append(island_list[i])
            i+=1
    return(final_islands)



# first identifies intersecting islands and then uses diff function to find UNintersected islands
def find_unintersected_islands(island_list, island_list_2):

    intersection_islands = []
    i = 0
    second_island_beginning = island_list_2[i][1]
    second_island_end = island_list_2[i][2]

    chrom_name_old_first_isl = island_list[0][0]
    # whole second island before first island flag
    flag = 0

    for (j,island) in enumerate(island_list):

        # if no more islands in control are left
        if i>=len(island_list_2):
            continue

        first_island_beginning = island[1]
        first_island_end = island[2]


        # check the chromosome change:
        # chrom names for chromosome switching
        chrom_name_current_first_isl = island[0]
        chrom_name_current_second_isl = island_list_2[i][0]

        if flag == 1:
            if chrom_name_current_first_isl != chrom_name_old_first_isl:
                chrom_name_old_first_isl = chrom_name_current_first_isl
                flag = 0
            else:
                continue
        # chromosome in first island list changes
        if chrom_name_current_first_isl != chrom_name_old_first_isl:
            chrom_name_old_first_isl = chrom_name_current_first_isl
            # go to next second island until chr names match
            while True:
                if ((chrom_name_current_second_isl == chrom_name_current_first_isl)) or (i>= len(island_list_2)-1):
                    break
                else:
                    i+=1
                    chrom_name_current_second_isl = island_list_2[i][0]

        while True:
            if i>=len(island_list_2):
                break

            # Chromosome change happens while we iterate through islands from the second list
            chrom_name_current_second_isl = island_list_2[i][0]
            if chrom_name_current_second_isl != chrom_name_current_first_isl:
                flag = 1
                break


            second_island_beginning = island_list_2[i][1]
            second_island_end = island_list_2[i][2]

            # Evaluating second list vs first list  island positions

            # Intersection condition: this one can be easily changed if needed
            if((second_island_beginning>=first_island_beginning) and (second_island_beginning<=first_island_end)) or\
                        ((second_island_end>=first_island_beginning) and (second_island_end<=first_island_end)) or \
                        ((second_island_beginning<first_island_beginning) and (second_island_end>first_island_end)):
                    intersection_islands.append(island)
                    i+=1
                    break
            # Intersection condition is not satisfied and whole second island is before whole first one: take next second island
            elif (second_island_end < first_island_beginning):
                i +=1
            # Intersection condition is not satisfied and whole second island is before whole first one: break and take next first island
            elif (second_island_beginning > first_island_end):
                break
    final_islands = diff(island_list,intersection_islands)

    return (final_islands)
"""