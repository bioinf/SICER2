import scipy, numpy, scipy.stats

# ---------------------------------------------------------------------------- #
print 'Step 1 of 3: Counting unique reads'

track = bgunzip(args.track.name)
reads, length, info, f_size = count_unique_reads(track, 1)

control = False
if args.control != None : control = bgunzip(args.control.name)

effective_length = count_effective_length(gms, track)
track.file_handle.close()

if args.window > max_window : args.window = max_window
if args.window < 10 : args.window = 10

lambdaa = count_lambda(reads, args.window, effective_length)

pvalue = args.pvalue
if pvalue > 1 or pvalue <= 0 : pvalue = 0.2
l0 = scipy.stats.poisson.ppf(1 - pvalue, lambdaa)

msg = "Window read threshold is {} reads, " + \
      "\ni.e. {} is minimum number of reads in window to consider" + \
      "\nthis window `eligible` with Poisson distribution p-value {}\n"
logging.info(msg.format(l0, l0, pvalue))

# ---------------------------------------------------------------------------- #

print 'Step 2 of 3: Making window list'

if args.window < 1 : args.window = 1
args.gap = args.gap/args.window
if args.gap < 0 : args.gap = 0

NORM_CONSTANT = 1
info = make_windows_list(info, l0, args.window, args.gap, NORM_CONSTANT)
logging.info('-' * 80)

# ---------------------------------------------------------------------------- #

print 'Step 3 of 3: Writing found islands'

if args.threshold < 0 : args.threshold = 0

# found, coverage = write_islands_list(info, lambdaa, args.window, l0, args.threshold, resultf)
found, coverage = write_islands(info, lambdaa, args.window, l0, args.threshold, resultf, f_size)
coverage_txt = str(coverage) + 'bp'
if coverage > 1000000 :
  coverage_txt = str(float(coverage)/1000000) + 'Mbp (' + coverage_txt + ')'
logging.info("Coverage: {}".format(coverage_txt))

msg = "There are {} islands found"
logging.info(msg.format(found))
print msg.format(found)

logging.info(('\n ').join([
 "Column names in BED file:",
 "1. Chromosome name",
 "2. Island start",
 "3. Island end",
 "4. Island score",
 "5. Number of reads in the island",
 "6. Number of eligible windows per island",
 "7. Number of gaps in the island"
]))

msg = "Finished. Elapsed time, minutes: {}"
logging.info(msg.format((time.time() - startTime) / 60))
logging.info("")


