from math import *
import scipy, numpy, scipy.stats

# ---------------------------------------------------------------------------- #
print '-' * 80
print 'Step 1 of 3: Counting unique reads'

track = bgunzip(args.track.name)
data, length, fragment, plambda, total_reads, effective_len = preparation(track, gms)
track.file_handle.close()

if args.window > max_window : args.window = max_window
if args.window < 10 : args.window = 10

pvalue = args.pvalue
if pvalue > 1 or pvalue <= 0 : pvalue = 0.2
l0 = scipy.stats.poisson.ppf(1 - pvalue, plambda)

msg = "Window read threshold:   {}\n" + \
      "i.e. {} is minimum number of reads in window to consider\n" + \
      "this window `eligible` with Poisson distribution p-value {}\n"
logging.info(msg.format(l0, l0, pvalue))

# ---------------------------------------------------------------------------- #

print 'Step 2 of 3: Making window list'

if args.window < 1 : args.window = 1
wlist = windows(data, length, fragment, args.window, gap_size)

# ---------------------------------------------------------------------------- #

print 'Step 3 of 3: Writing found islands'

if threshold > 0 :
  logging.info("The score threshold is: {}".format(threshold))
else :
  e_value = 100.0; bin_size = 0.001
  bg = Island_threshold(total_reads, args.window, gap_size, pvalue, effective_len, bin_size, e_value)
  threshold = bg.threshold
  logging.info("The score threshold is: {} (auto define)".format(threshold))

found, coverage = islands(wlist, l0, plambda, args.window, threshold, resultf)
coverage_txt = str(coverage) + ' bp'
if coverage > 1000000 :
  coverage_txt = str(float(coverage)/1000000) + ' Mbp (' + coverage_txt + ')'
logging.info("Coverage: {}".format(coverage_txt))
logging.info("Islands found: {}".format(found))

msg = "Done.\nIslands found:    {}\nOutput file name: {}\nLog file name:    {}"
print msg.format(found, resultf, logfile)

logging.info(('\n ').join([
  "Output file:",
  "1. Chromosome name",
  "2. Island start",
  "3. Island end",
  "4. Island name",
  "4. Island score"
]))

msg = "\nFinished. Elapsed time, seconds: {}"
logging.info(msg.format(round(time.time() - startTime, 2)))
logging.info("\n")
