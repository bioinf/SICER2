#!/usr/bin/env python
# -*- coding: utf-8 -*-

import scipy, scipy.stats, numpy, pysam, os.path, argparse, logging

def main():
  parser = argparse.ArgumentParser(
    prog        = "SICER2", 
    description = "Tool for ChiP-seq analysis to find broad peaks"
  )
  print "Ok"


def get_chromosomes_info(bam):
  # Есть ли индекс бам-файла
  bai = bam + ".bai"
  if not os.path.isfile(bai):
    # No BAM index file was found, new index was generated
    pysam.index(bam)
    print "# No BAM index file was found, new index was generated"
  # Общая информация
  print "Collecting information about sample from .bai file:"
  print "[ref.seq. name, ref.seq. length, number of mapped and unmapped reads]"
  chromosomes_info = []
  try:
    for chr in pysam.idxstats(bam):
      chromosomes_info.append(chr.split("\t")[:-1])
    chromosomes_info.pop() # Last line is unmapped reads, we don't need them
  except:
    print "PROBLEM WITH BAM FILE OR pysam.idxstats() COMMAND\nYour BAM file {} probably is not sorted."
    print "To sort it with samtools use comand: \n'samtools sort ... "
  return chromosomes_info

def count_unique_reads(bam_path, chromosomes_info):
  bamfile = pysam.AlignmentFile(bam_path, 'rb')
  total_reads_count = 0
  total_unique_reads_count = 0
  previous_read_strand = 0
  plus_reads_count = 0
  minus_reads_count = 0
  all_read_length = set()
  for chromosome in chromosomes_info:
    chr_unique_reads_count = 0
    chr_total_reads_count = 0
    beginning_of_the_previous_read = 0
    current_chromosome_name = chromosome[0]
    # currentChromosomeSize = int(chromosome[1])
    all_reads_in_chromosome = bamfile.fetch(current_chromosome_name)
    for read in all_reads_in_chromosome:
      read_str = str(read)
      # read strand: 0 = +         16 = -
      read_strand = ([int(s) for s in read_str.split() if s.isdigit()][0])
      beginning_of_the_read = ([int(s) for s in read_str.split() if s.isdigit()][2])
      # print read_str
      # print read_strand
      if beginning_of_the_read != beginning_of_the_previous_read or read_strand != previous_read_strand:
        beginning_of_the_previous_read = beginning_of_the_read
        if read_strand == 0:
          minus_reads_count += 1
        else:
          plus_reads_count += 1
        previous_read_strand = read_strand
        all_read_length.add(len(read_str.split()[9]))
        # print read_str.split()[9]
        # print ''
        total_unique_reads_count += 1
        chr_unique_reads_count += 1
      chr_total_reads_count += 1
    print("On {} there are {} unique reads among {}".format(current_chromosome_name,
                                                                   chr_unique_reads_count,
                                                                   chr_total_reads_count,))
    total_reads_count += chr_total_reads_count
  # print("Unique reads counted")
  bamfile.close()
  if len(all_read_length) == 1:
      print("\nAverage read length {} bp".format(all_read_length.pop()))
  else:
      print('\nVariable read length {} bp'.format([_ for _ in all_read_length]))
  print("Library depth: there are {} unique reads out of {}.\nIn other words {} % of reads are unique".
               format(total_unique_reads_count, total_reads_count,
                      round(float(total_unique_reads_count)/float(total_reads_count)*100, 1)))
  print("Strand symmetry: \n {} (+) \n {} (-)".format(plus_reads_count, minus_reads_count))
  return total_unique_reads_count

# i = get_chromosomes_info('data/chrM.bam')
# count_unique_reads('data/chrM.bam', i)


def count_effective_length(effective_proportion, chromosomes_info):
    total_genome_length = sum(int(row[1]) for row in chromosomes_info)
    effective_length = effective_proportion * total_genome_length
    return effective_length

def count_lambda(unique_reads_count, window_size, effective_length):
    lambdaa = float(window_size) * float(unique_reads_count) / float(effective_length)
    print("\nAverage density of reads per {} bp window is {}".format(window_size, round(lambdaa, 2)))
    return lambdaa

def make_windows_list(bam_path, chromosomes_info, l0, window_size, gap, unique_reads_count, normalization_coef):
    print("Making eligible windows of {} bp with allowed gap_size {} bp".format(window_size, window_size*gap))
    bamfile = pysam.AlignmentFile(bam_path, 'rb')
    window_list = []
    window_list_dict = {}
    # print("chromosome_name, chromosome_size, total_number_of_eligible_windows_on_chromosome")
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
        print("On {} there are {} eligible windows".format(current_chromosome_name, chr_window))
    # candidate windows == eligible + ineligible(below gap)
    print("\nThere are {} candidate windows".format(len(window_list) - i))
    window_list.append([1, 1])
    bamfile.close()
    # print(window_list)
    # print(window_list_dict)
    return (window_list, window_list_dict)


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
    print("There are {} islands found".format(len(islands_list)))
    #print(zero_score_islands)
    # print(islands_list)
    return islands_list


def broadpeaks_wo_control(bam_path, WINDOW_SIZE, GAP, EFFECTIVE_PROPORTION, ISLAND_SCORE_THRESHOLD, p0):
    chromosomes_info = get_chromosomes_info(bam_path)

    print("\nStep 1 of 4\nCOUNTING UNIQUE READS\n")
    input_unique_reads_count = count_unique_reads(bam_path, chromosomes_info)

    # Effective genome length (L)
    effective_length = count_effective_length(EFFECTIVE_PROPORTION, chromosomes_info)

    # Lambda for poisson distribution
    lambdaa = count_lambda(input_unique_reads_count, WINDOW_SIZE, effective_length)

    # Minimum #reads in a window for eligibility
    # Formula (1), finding l0
    l0 = scipy.stats.poisson.ppf(1 - p0, lambdaa)
    NORMALIZATION_CONSTANT = 1
    print("\nWindow read threshold is {} reads, \ni.e. {} is minimum number of reads in window "
                 "to consider this window `eligible` with Poisson distribution p-value {}".format(l0, l0, p0))

    print("\nStep 2 of 4\nMAKING WINDOW LIST\n")

    (window_list, window_list_dict) = make_windows_list(bam_path, chromosomes_info, l0, WINDOW_SIZE, GAP,
                                            input_unique_reads_count, NORMALIZATION_CONSTANT)
    print("\nStep 3 of 4\nMAKING ISLAND LIST\n")

    island_list = make_islands_list(window_list, lambdaa, WINDOW_SIZE, l0, chromosomes_info,
                                            ISLAND_SCORE_THRESHOLD)

    return island_list


island_list = broadpeaks_wo_control('data/H3K36me3_GSM1261677.bam', 200, 1, 0.77, 0, 0.01)

print("\nStep 4 of 4\nWRITING FOUND ISLANDS TO `{}` BED FILE\n".format('src'))

print len(island_list)
print("Done")


