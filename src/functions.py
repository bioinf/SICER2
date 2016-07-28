# ---------------------------------------------------------------------------- #
def define_specie(info) :
  vote = {}
  total = 1
  for c in info :
    total += info[c]['Unique reads']

  for c in info :
    dists = []
    for table in species :
      for chr in species[table]:
        dists.append([pow(chr - info[c]['Length'], 2), table])
    v = sorted(dists, key=lambda e: e[0])
    for k in [0,1,2] :
      if v[k][1] not in vote : vote[v[k][1]] = 0
      vote[v[k][1]] += float(3 - k) * info[c]['Unique reads']/total
  code = max(vote, key = vote.get)
  return [code, spnames[code]]


# ---------------------------------------------------------------------------- #
def parser(data_generator):
  chunk = next(data_generator)
  CtoPy = { 'A':'<c', 'c':'<b', 'C':'<B', 's':'<h', 'S':'<H', 'i':'<i', 'I':'<I' }
  py4py = { 'A':  1,  'c':  1,  'C':  1,  's':  2,  'S':  2,  'i':  4 , 'I':  4  }
  dna = '=ACMGRSVTWYHKDBN'
  cigar_codes = 'MIDNSHP=X'
  from array import array
  from struct import unpack
  p = 0
  while True:
    try:
      while len(chunk) < p+36: chunk = chunk[p:] + next(data_generator); p = 0
      block_size,refID,pos,l_read_name,mapq,bin_,n_cigar_op,flag,l_seq                          = unpack('<iiiBBHHHi'   ,chunk[p:p+24])
      while len(chunk) < p + 4 + block_size: chunk = chunk[p:] + next(data_generator); p = 0
    except StopIteration: break
    end = p + block_size + 4
    p = end
    yield pos,refID,flag,l_seq


# ---------------------------------------------------------------------------- #
def chrsort(s) :
  s = s.lower()  
  if s[0:3] == 'chr' : 
    s = s.replace('chr','')
    if unicode(s, 'utf-8').isnumeric() : return int(s)
    return 40
  return 999

def count_unique_reads(track, maxdr) :
  reads = { 0 : 0, 16 : 0 }
  drate = { 0 : 0, 16 : 0 }
  lhist = {}
  info = {}
  
  for read in parser(track) :
    pos, chr, strand, l_seq = read
    if chr < 0 or chr > len(track.chromosome_names) : continue
    if l_seq not in lhist : lhist[l_seq] = 0
    lhist[l_seq] += 1
    c = track.chromosome_names[chr]
    # New chromosome name:
    if c not in info : 
      info[c] = {
        'Data' : array.array('l', []), # position * 10 + strand indicator (0 or 1)
        'Length' : track.chromosome_lengths[chr],
        'Unique reads' : 0,
        'Total reads' : 0,
        'Beginning of the previous' : 0
      }
    # Два рида в одном месте
    if pos == info[c]['Beginning of the previous'] :
      drate[strand] += 1
      if drate[strand] <= maxdr :
        info[c]['Data'].append(10 * pos + (strand % 15))
        info[c]['Unique reads'] += 1
        reads[strand] += 1
    else :
      drate = { 0 : 0, 16 : 0 }
      info[c]['Beginning of the previous'] = pos
      info[c]['Data'].append(10 * pos + (strand % 15))
      info[c]['Unique reads'] += 1
      reads[strand] += 1
    info[c]['Total reads'] += 1
  # <- for

  total_reads = 0; unique_reads = 0
  logging.info("Chromosome, unique reads, total reads:")

  keys = info.keys()
  keys = sorted(keys, key = lambda (c): chrsort(c))
  for c in keys :
    total_reads  += info[c]['Total reads']
    unique_reads += info[c]['Unique reads']
    logging.info("{}\t{}\t{}".format(c, info[c]['Unique reads'], info[c]['Total reads']))
  # <- for
  
  mean_length = 0
  for l in lhist :
    mean_length += l * lhist[l]
  mean_length = mean_length/total_reads

  logging.info('-' * 80)
  
  msg = "Average read length: {}"
  logging.info(msg.format(mean_length))

  msg = "Library depth: there are {} unique reads out of {}"
  logging.info(msg.format(unique_reads, total_reads))

  msg = "Duplication rate: {}"
  logging.info(msg.format(round(1 - float(unique_reads)/float(total_reads)*100, 1)))

  msg = "Strand symmetry:\n  {} (+)\n  {} (-)"
  logging.info(msg.format(reads[0], reads[16]))

  return [unique_reads, mean_length, info]


# ---------------------------------------------------------------------------- #
def get_gms(info, length, gms) :
  if gms > 0 and gms <= 1 :
    msg = "Using custom effective proportion: {}"
    logging.info(msg.format(gms))
    return gms
  
  specie, specie_name = define_specie(info)
  GMS = {}
  GMS['hg19'] = [[20,0.697969],[25,0.750921],[30,0.782328],[40,0.824712],[60,0.867066],[80,0.882649],[100,0.889473],[120,0.893435],[140,0.896156]]
  GMS['mm10'] = [[25,0.766128],[30,0.790964],[40,0.822389],[60,0.857177],[80,0.875752],[100,0.88708],[120,0.894755],[140,0.900449]]

  # Default value
  if specie not in GMS : 
    gms = 0.77
    msg = "Using default effective proportion: {}"
    logging.info(msg.format(gms))
    return gms

  # Approximation
  if GMS[specie][0][0] > length :
    gms = GMS[specie][0][1]
  if GMS[specie][-1][0] < length :
    gms = GMS[specie][-1][1]
  
  for e in range(0, len(GMS[specie]) - 1) :
    a1, v1 = GMS[specie][e]
    a2, v2 = GMS[specie][e + 1]
    if a2 >= length and a1 <= length :
      # (a2 - a1)/(v2 - v1) == (length - a1)/(x - v1)
      gms = (length - a1) * (v2 - v1)/(a2 - a1) + v1

  # ...
  msg = "Using default effective proportion for {}: {}"
  logging.info(msg.format(specie_name, gms))
  return gms


# ---------------------------------------------------------------------------- #
def count_effective_length(effective_proportion, track) :
  return effective_proportion * sum(track.chromosome_lengths)


# ---------------------------------------------------------------------------- #
def count_lambda(unique_reads_count, wsize, effective_length):
  lambdaa = float(wsize) * float(unique_reads_count) / float(effective_length)
  msg = "Average density of reads per {} bp window is {}"
  logging.info(msg.format(wsize, round(lambdaa, 2)))
  return lambdaa


# ---------------------------------------------------------------------------- #
def make_windows_list(info, l0, window_size, gap, normalization_coef):
  msg = "Making eligible windows of {} bp with allowed gap_size {} bp"
  logging.info(msg.format(window_size, window_size * gap))
  logging.info('-' * 80)
  logging.info("Chromosome, eligible windows:")

  keys = info.keys()
  keys = sorted(keys, key = lambda (c): chrsort(c))
  for c in keys :
    info[c]['Windows'] = array.array('l', [])
    previous = 0
    previous_read_strand = 0
    gap_count = 0
    window_start = 0
    window_reads_count = 0
    chr_window = 0
    for read in info[c]['Data'] :
      read_strand = read % 10
      begin = (read - read_strand) / 10
      if (begin != previous) or (read_strand != previous_read_strand):
        previous = begin
        previous_read_strand = read_strand
        gap_flag = True
        while True:
          if window_start <= begin < window_start + window_size:
            window_reads_count += 1
            break
          elif begin < window_start:
            break
          else:
            window_reads_count = int(float(window_reads_count) * normalization_coef)
            if window_reads_count < l0:
              gap_count += 1
            else:
              gap_flag = False
              gap_count = 0
            info[c]['Windows'].append(max_window * window_start + window_reads_count)
            chr_window += 1
            if gap_count > gap or gap_flag:
              gap_flag = True
              while gap_count > 0:
                info[c]['Windows'].pop()
                gap_count -= 1
                chr_window -= 1
            window_start += window_size
            window_reads_count = 0
    if window_reads_count != 0:
      info[c]['Windows'].append(max_window * window_start + window_reads_count)
    logging.info("{}\t{}".format(c, chr_window))
  return info

def write_islands_list(info, lambdaa, wsize, l0, threshold, resultf):
  f = open(resultf, 'w+'); islands = 0; coverage = 0

  def score(reads):
    if reads >= l0:
      temp = scipy.stats.poisson.pmf(reads, lambdaa)
      if temp < 1e-320: window_score = 1000
      else: window_score = -numpy.log(temp)
    else:
      window_score = 0
    return window_score

  def write(c, e) :
    if e['score'] < threshold : return 0, 0
    f.write(c+'\t'+str(e['from'])+'\t'+str(e['to'])+'\t'+str(e['score'])+'\t'+str(e['reads'])+'\t'+str(e['count'] - e['gaps'])+'\t'+str(e['gaps'])+'\n')
    return 1, (e['to'] - e['from'])

  def new_island(init, reads):
    return { 
      'from'  : init, 
      'to'    : init + wsize, 
      'score' : score(reads), 
      'reads' : reads, 
      'count' : 1, 
      'gaps'  : 0
    }

  for c in info :
    island = False
    for w in info[c]['Windows'] :
      reads = w % max_window
      init = (w - reads)/max_window
      if not island :
        island = new_island(init, reads)
      else :
        # Same island
        if init == island['to'] :
          island['to'] += wsize
          island['count'] += 1
          island['reads'] += reads
          island['score'] += score(reads)
          if reads == 0 : island['gaps'] += 1
        # Other island
        else :
          i, cov = write(c, island)
          islands += i; coverage += cov
          island = new_island(init, reads)
    # <- end for
    if island :
      i, cov = write(c, island)
      islands += i; coverage += cov
      # islands += write(c, island)

  # <- end for
  f.close()
  return islands, coverage
