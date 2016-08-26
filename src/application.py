import scipy, numpy, scipy.stats

# ---------------------------------------------------------------------------- #
print 'Step 1 of 3: Counting unique reads'

track = bgunzip(args.track.name)
data, length, fragment, plambda = preparation(track, gms)
track.file_handle.close()

if args.window > max_window : args.window = max_window
if args.window < 10 : args.window = 10

pvalue = args.pvalue
if pvalue > 1 or pvalue <= 0 : pvalue = 0.2
l0 = scipy.stats.poisson.ppf(1 - pvalue, plambda)

msg = "Window read threshold is {} reads, " + \
      "\ni.e. {} is minimum number of reads in window to consider" + \
      "\nthis window `eligible` with Poisson distribution p-value {}\n"
logging.info(msg.format(l0, l0, pvalue))

# ---------------------------------------------------------------------------- #

print 'Step 2 of 3: Making window list'

if args.window < 1 : args.window = 1
if args.gap < 1 : args.gap = 1

wlist = windows(data, length, args.window, gap_size)

# ---------------------------------------------------------------------------- #

print 'Step 3 of 3: Writing found islands'

if args.threshold < 0 : args.threshold = 0
found, coverage = islands(wlist, l0, plambda, args.window, args.threshold, resultf)

coverage_txt = str(coverage) + ' bp'
if coverage > 1000000 :
  coverage_txt = str(float(coverage)/1000000) + ' Mbp (' + coverage_txt + ')'
logging.info("Coverage: {}".format(coverage_txt))

logging.info("Islands found: {}".format(found))

logging.info(('\n ').join([
  "Output file:",
  "1. Chromosome name",
  "2. Island start",
  "3. Island end",
  "4. Length of island region",
  "5. * Absolute peak summit position",
  "6. * Pileup height at peak summit, -log10(pvalue) for the peak summit",
  "7. Island score",
  "8. Number of reads in the island",
  "9. Number of eligible windows per island",
  "10. Number of gaps in the island"
]))

msg = "Finished. Elapsed time, minutes: {}"
logging.info(msg.format((time.time() - startTime) / 60))
logging.info("")

