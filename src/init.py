#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys, array, argparse, logging, time, datetime

startTime = time.time()
app = os.path.basename(__file__)
max_window = 1000000

# Beautiful --help
class AParser(argparse.ArgumentParser) :
  def help_message(name = None) :
    print ''
    print 'Usage:'
    print '  ./' + app + ' -t [input bam-file] [Optional arguments]'
    print ''
    print 'Optional arguments:'
    print '  -t,  --input      Path to `input.bam` file'
    print '  -o,  --output     Output file name'
    print '  -l,  --log        Output log file name'
    print '  -w,  --window     Window size (in bp). Default: 200'
    print '  -f,  --fragment   Fragment size (in bp). Default: 250'
    print '  -e,  --gms        Proportion of effective genome length; has to be in [0.0, 1.0]. Default: auto'
    print '  -g,  --gap        Gap size shows how many bases could be skipped. Default: 200'
    print '  -p,  --pvalue     P-value; has to be in [0.0, 1.0]. Default: 0.2'
    print '  -x,  --threshold  Island score threshold. Default: auto'
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
parser.add_argument('-l', '--log',       default = None)
parser.add_argument('-c', '--control',   default = None, type = argparse.FileType('r'))
parser.add_argument('-w', '--window',    default = 200, type = int)
parser.add_argument('-f', '--fragment',  default = 'auto')
parser.add_argument('-e', '--gms',       default = 'auto')
parser.add_argument('-g', '--gap',       default = 200, type = int)
parser.add_argument('-p', '--pvalue',    default = 0.2, type = float)
parser.add_argument('-x', '--threshold', default = 'auto')

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

if args.log : logfile = args.log
if args.output : resultf = args.output.name #  + '_peaks.bed'

logging.basicConfig(filename = logfile, level = logging.DEBUG, format = '%(message)s')
logging.info('=' * 80)
logging.info('Date: ' + datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'))
logging.info('Command:')

cmd = ['--track ' + str(args.track.name)]

try :
  fragmentsize = float(args.fragment)
except :
  fragmentsize = 0

try :
  threshold = float(args.threshold)
except :
  threshold = 0

if args.control : cmd.append('--control ' + str(args.control.name))
if args.output  : cmd.append('--output ' + str(args.output.name))
if args.log     : cmd.append('--log ' + str(args.log))
if fragmentsize > 0 : cmd.append('--fragment ' + str(fragmentsize))

try :
  gms = float(args.gms)
except :
  gms = 0

gap_size = max(args.gap, args.window)
cmd.extend([
  '--window ' + str(args.window),
  '--gms ' + str(args.gms),
  '--gap ' + str(gap_size),
  '--pvalue ' + str(args.pvalue),
  '--threshold ' + str(args.threshold)
])

logging.info(sys.argv[0] + ' \\\n  ' + (' \\\n  ').join(cmd))
logging.info('')
