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
    p += 36
    qname = chunk[p:p+l_read_name-1]
    p = end
    yield pos,refID,flag,l_seq,qname


# ---------------------------------------------------------------------------- #
def chrsort(s) :
  s = s.lower()  
  if s[0:3] == 'chr' : 
    s = s.replace('chr','')
    if unicode(s, 'utf-8').isnumeric() : return int(s)
    return 40
  return 999


# ---------------------------------------------------------------------------- #
def beautiful_table(data, cols = 3) :
  t = lambda A : [list(x) for x in zip(*A)] # Транспонирование
  r = (len(data) - len(data)%cols)/cols + 1 # Число строк в столбце
  T = []
  for i in range(0, len(data), r):
    col = data[i:i + r]
    while len(col) <= r : 
      col.append(len(col[0]) * [' '])
    col = t(col)
    siz = map(lambda c : max(map(lambda x : len(str(x)), c)), col)
    # ¯\_(ツ)_/¯
    for i in range(len(col)) :
      k = map(lambda x : str(x) + (siz[i] - len(str(x)) + 1) * " ", col[i])
      T.append(k)
    T.append(['  |  '] * r)
  if len(T) == 0 : return ''
  T = map(lambda x : ('').join(x), t(T[0:-1]))
  V = (len(T[0])) * '-' + '\n'
  return V + ('\n').join(T) + '\n' + V


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
      # Интерполяция
      # (a2 - a1)/(v2 - v1) == (length - a1)/(x - v1)
      gms = (length - a1) * (v2 - v1)/(a2 - a1) + v1

  # ...
  msg = "Using default effective proportion for {}: {}"
  logging.info(msg.format(specie_name, gms))
  return gms


# ---------------------------------------------------------------------------- #
def count_lambda(unique_reads_count, wsize, effective_length):
  plambda = float(wsize) * float(unique_reads_count) / float(effective_length)
  msg = "Average density of reads per {} bp window: {}"
  logging.info(msg.format(wsize, round(plambda, 2)))
  return plambda


# ---------------------------------------------------------------------------- #
def preparation(track, gms) :
  length_hist = {}
  distance = { 'sum' : 0, 'count' : 0 }
  reads = { 0 : 0, 16 : 0 }
  info = {}
  data = {}

  for read in parser(track) :
    pos, chr, strand_key, l_seq, qname = read
    if chr < 0 or chr > len(track.chromosome_names) : continue
    c = track.chromosome_names[chr]
    if l_seq not in length_hist : length_hist[l_seq] = 0
    length_hist[l_seq] += 1
    # pos += l_seq/2

    # New chromosome name:
    if c not in info : 
      firsts = 0
      data[c] = array.array('l', [])
      info[c] = {
        'Length' : track.chromosome_lengths[chr],
        'Unique reads' : 0,
        'Total reads'  : 0,
        'Names' : {}
      }

    firsts += 1
    if firsts < 300 :
      if strand_key in [99,147,163,83] : # check paired reads
        if qname not in info[c]['Names'] : 
          info[c]['Names'][qname] = pos
        else :
          dist = abs(info[c]['Names'][qname] - pos)
          if dist > 0 :
            distance['sum'] += dist
            distance['count'] += 1

    info[c]['Total reads'] += 1
    if strand_key not in reads : reads[strand_key] = 0
    reads[strand_key] += 1

    prev_key = 'Previous ' + str(strand_key)
    ignore = False
    if prev_key in info[c] :
      if info[c][prev_key] == pos : ignore = True

    if ignore == False :
      info[c][prev_key] = pos
      data[c].append(1000 * pos + strand_key)
      info[c]['Unique reads'] += 1
    # <- if
  # <- for read in parser(track)

  total_reads = 0; unique_reads = 0
  keys = info.keys()
  keys = sorted(keys, key = lambda (c): chrsort(c))
  
  tbl = []
  for c in keys :
    tbl.append([c, info[c]['Unique reads'], info[c]['Total reads']])
    total_reads  += info[c]['Total reads']
    unique_reads += info[c]['Unique reads']
  logging.info("Chromosome name, Unique reads, Total reads:")
  logging.info(beautiful_table(tbl))

  mean_length = 0
  for l in length_hist :
    mean_length += l * length_hist[l]
  mean_length = mean_length/total_reads

  msg = "Average read length: {}"
  logging.info(msg.format(mean_length))

  drate = round((1 - float(unique_reads)/float(total_reads)) * 100, 1)
  msg = "\nLibrary depth\n  Duplication rate:  {}%\n  Total reads:       {}\n  Unique reads:      {}"
  logging.info(msg.format(drate, total_reads, unique_reads))
  
  notes = [
    [99 , "  - Reads mapped in proper pair. Mate reverse strand, first in pair:  {}"], # +
    [147, "  - Reads mapped in proper pair. Read reverse strand, second in pair: {}"], # +
    [83 , "  - Reads mapped in proper pair. Read reverse strand, first in pair:  {}"], # -
    [163, "  - Reads mapped in proper pair. Mate reverse strand, second in pair: {}"], # -

    [97 , "  - Mate reverse strand, first in pair:                               {}"], #
    [161, "  - Mate reverse strand, second in pair:                              {}"], # -
    [81 , "  - Read reverse strand, first in pair:                               {}"], # +
    [145, "  - Read reverse strand, second in pair:                              {}"], #
  			 
    [113, "  - Mate reverse strand, Read reverse strand, first in pair:          {}"],
    [177, "  - Mate reverse strand, Read reverse strand, second in pair:         {}"],
    [65 , "  - First in pair:                                                    {}"],
    [129, "  - Second in pair:                                                   {}"],
  			 
    [73 , "  - Mate unmapped, first in pair:                                     {}"],
    [137, "  - Mate unmapped, second in pair:                                    {}"],
    [89 , "  - Mate unmapped, read reverse strand, first in pair:                {}"],
    [153, "  - Mate unmapped, read reverse strand, second in pair:               {}"]
  ]

  def paired() :
    for x in notes :
      if x[0] in reads : return True
    return False

  logging.info("\nStrand symmetry")
  fragment_size = fragmentsize
  
  if paired() :
    if fragment_size == 0 : 
      fragment_size = distance['sum']/distance['count'] + mean_length
    for key, msg in notes :
      logging.info(msg.format(reads[key] if key in reads else 0))
  else :
    if fragment_size == 0 : 
      fragment_size = 250
    logging.info("  [+] " + str(reads[0]))
    logging.info("  [-] " + str(reads[16]))

  logging.info("\nFragment size: " + str(fragment_size))
  logging.info("")
  
  effective_len = get_gms(info, mean_length, gms) * sum(track.chromosome_lengths)
  plambda = count_lambda(unique_reads, args.window, effective_len)

  logging.info("")
  logging.info("Genome Length:           {}".format(sum(track.chromosome_lengths)))
  logging.info("Effective genome Length: {}".format(effective_len))

  return [data, mean_length, fragment_size, plambda, total_reads, effective_len]


# ---------------------------------------------------------------------------- #
def windows(data, length, fragment, window_size, gap_size, l0):
  keys = sorted(data.keys(), key = lambda (c): chrsort(c))
  tbl = []; eligible = 0; gaps = [0,0,0,0] # gaps - histogram
  for c in keys :
    # wlist = array.array('l', [])
    wlist = []
    last_init = -max_window
    chr_windows = 0
    for read in data[c] :
      strand, init = read%1000, read/1000
      gp = (init - last_init)/window_size
      if gp <= 3 : gaps[gp] += 1

      w_init = init/window_size * window_size
      w_last_init = last_init/window_size * window_size

      if w_last_init == w_init :
        wlist[-1] += 1
        chr_windows += 1
      else :
        wlist.append(w_init * max_window + 1)
        chr_windows += 1
      last_init = init
    # <- for read in data[c]

    wlist_good = [] # array.array('l', [])
    for i in wlist :
      if i%max_window >= l0 :
        wlist_good.append(i)

    last_init = -max_window
    wlist = array.array('l', [0])
    for i in wlist_good :
      reads, init = i%max_window, i/max_window
      if last_init + window_size + gap_size > init:
         gg = (init - last_init - window_size)/window_size
         for k in range(gg) :
           wlist.append((last_init + (k + 1) * window_size) * max_window)
      wlist.append(i)
      last_init = init
    # <- FOR

    tbl.append([c, chr_windows])
    eligible += chr_windows
    data[c] = wlist
  msg = "Total eligible windows of {}bp with allowed gap size {}bp: {}"
  logging.info(msg.format(window_size, gap_size, chr_windows))

  msg = "Chromosome name, Eligible windows:\n{}"
  logging.info(msg.format(beautiful_table(tbl)))

  msg = "Gap size count:"
  for i in range(4) :
    msg += "\n  {:>3} - {:>3}bp: {}".format(i * window_size + 1, (i + 1) * window_size, gaps[i])
  logging.info(msg)

  return data

# ---------------------------------------------------------------------------- #
def islands(wlist, l0, plambda, window_size, threshold, resultf, fragment):
  f = open(resultf, 'w+')
  islands = 0; coverage = 0; reads_scores = {}

  def score(reads):
    if reads <= l0 : return 0
    t = scipy.stats.poisson.pmf(reads, plambda)
    return 1000 if t < 1e-320 else -numpy.log(t)

  def new_island(init, reads):
    return { 
      'from'  : init,
      'to'    : init + window_size,
      'score' : reads_scores[reads],
      'reads' : reads, 
      'count' : 1, 
      'gaps'  : 0
    }

  def write(chromosome, e, islands) :
    if e['score'] < threshold: return 0, 0
    island_size = e['to'] - e['from']#  - fragment
    if island_size < 30 : return 0, 0
    name = 'is_' + str(islands)
    # f.write(chromosome +'\t' + str(e['from'] + fragment/2) +'\t' + str(e['to'] - fragment/2)+'\t' + name +'\t' + str(e['score']) + '\n')
    f.write(chromosome +'\t' + str(e['from']) +'\t' + str(e['to'])+'\t' + name +'\t' + str(e['score']) + '\n')
    return 1, island_size

  for c in wlist :
    island = False
    for w in wlist[c] :
      reads, init = w%max_window, w/max_window
      # if reads < l0 : continue
      if reads not in reads_scores : 
        reads_scores[reads] = score(reads)
      if not island :
        island = new_island(init, reads)
      else :
        # Same island
        if init == island['to'] :
          island['to'] += window_size
          island['count'] += 1
          island['reads'] += reads
          island['score'] += reads_scores[reads]
          if reads == 0 : island['gaps'] += 1
        # Other island
        else :
          i, cov = write(c, island, islands)
          islands += i; coverage += cov
          island = new_island(init, reads)
    # <- end for
    if island :
      i, cov = write(c, island, islands)
      islands += i; coverage += cov
  # <- end for
  f.close()
  return islands, coverage
