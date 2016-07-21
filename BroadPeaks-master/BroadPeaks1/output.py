#!/usr/bin/env python
import logging


def write_output(outfile, island_list, island_score_threshold):
    f = open(outfile, 'wb')
    logging.info("Line in BED file: "
                 "'chromosome_name\tisland_start\tisland_end\tisland_score\tnumber_of_reads_in_the_island\t"
                 "number_of_eligible_windows_per_island\tnumber_of_gaps_in_the_island\tscore_for_fdr\tfdr'")
    for island in island_list:
        if island[3] >= island_score_threshold:
            island_string = '\t'.join([str(_) for _ in island]) + '\n'
# str(island[0]) + "\t" + str(island[1]) + "\t" + str(island[2]) + "\t" + str(island[3]) + "\t" + str(island[4]) + "\t" + str(island[5])+ "\t" + str(island[6]) + "\n"
            f.write(island_string)
    # print(island_list[11104])
    f.close()
