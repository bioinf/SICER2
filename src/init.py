#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, array, argparse, logging, time, datetime

startTime = time.time()
app = os.path.basename(__file__)
max_window = 1000000

# Красивый хелп
class AParser(argparse.ArgumentParser) :
  def help_message(name = None) :
    print ''
    print 'Usage:'
    print '  ./' + app + ' -t [input bam-file] [Optional arguments]'
    print ''
    print 'Optional arguments:'
    print '  -t,  --input      Path to `input.bam` file'
    print '  -o,  --output     Output file name'
    # print '  -c,  --control    Path to `control.bam` file.  DEFAULT: no control file'
    print '  -w,  --window     Window size (bp).  DEFAULT: 200'
    print '  -e,  --gms        Proportion of effective genome length; has to be in [0.0, 1.0]  DEFAULT: auto'
    print '  -g,  --gap        Gap size shows how many bases could be skipped  DEFAULT: 200'
    print '  -p,  --pvalue     P-value; has to be in [0.0, 1.0]  DEFAULT: 0.1'
    print '  -x,  --threshold  Island score threshold  DEFAULT: 0'
    print ''
    print 'Examples:'
    print '  ./' + app + ' input.bam'
    print '  ./' + app + ' -t input.bam -c control.bam -w 100 -g 100 -e 0.845'
    print ''
  def error(self, message):
    sys.stderr.write('Error:\n  %s\n' % message)
    self.help_message()
    sys.exit(2)

# ---------------------------------------------------------------------------- #

parser = AParser()
parser.add_argument('-t', '--track',     type = argparse.FileType('r'))
parser.add_argument('-o', '--output',    type = argparse.FileType('w'))
parser.add_argument('-c', '--control',   default = None, type = argparse.FileType('r'))
parser.add_argument('-w', '--window',    default = 200, type = int)
parser.add_argument('-e', '--gms',       default = 'auto')
parser.add_argument('-g', '--gap',       default = 200, type = int)
parser.add_argument('-p', '--pvalue',    default = 0.1, type = float)
parser.add_argument('-x', '--threshold', default = 0, type = float)

if len(sys.argv) == 1 :
  parser.help_message()
  sys.exit(2)

if len(sys.argv) == 2 :
  args = parser.parse_args(['-t', sys.argv[1]])
else :
  args = parser.parse_args()

if args.track.name[-4:] == '.bam' :
  logfile = args.track.name[0:-4] + '_output.log'
  resultf = args.track.name[0:-4] + '_peaks.bed'
else :
  logfile = args.track.name + '_output.log'
  resultf = args.track.name + '_peaks.bed'
if args.output :
  resultf = args.output.name #  + '_peaks.bed'

logging.basicConfig(filename = logfile, level = logging.DEBUG, format = '%(message)s')
logging.info('=' * 80)
logging.info('Date: ' + datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
logging.info('Command:')

cmd = ['--track ' + str(args.track.name)]
if args.control : cmd.append('--control ' + str(args.control.name))
if args.output : cmd.append('--output ' + str(args.output.name))

cmd.extend([
  '--window ' + str(args.window),
  '--gms ' + str(args.gms),
  '--gap ' + str(args.gap),
  '--pvalue ' + str(args.pvalue),
  '--threshold ' + str(args.threshold)
])
logging.info(sys.argv[0] + ' \\\n  ' + (' \\\n  ').join(cmd))
logging.info('-' * 80)
