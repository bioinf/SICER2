# Source:
# https://github.com/JohnLonginotto/pybam

import zlib
import struct
import itertools
import subprocess

class bgunzip:
    def __init__(self,file_handle,blocks_at_a_time=50):
        ## First we init some basic things so that we will fill up in a sec.
        self.blocks_at_a_time = blocks_at_a_time
        self.bytes_read = 0
        self.rows = 0
        self.chromosome_names = []
        self.chromosome_lengths = []
        self.chromosomes_from_header = []
        self.file_handle = file_handle

        # Now we init the generator that, when iterated, will return some chunk of uncompress data.
        self.generator = self._generator()

        # We grab the first chunk, which is more than enough to contain all the header information, and parse the binary header:
        first_chunk = next(self.generator)
        if first_chunk[:4] != 'BAM\1': print 'ERROR: This doesnt look like a bam file to me :/'; exit()
        length_of_header = struct.unpack('<i',first_chunk[4:8])[0]
        self.header_text = struct.unpack(str(length_of_header)+'s',first_chunk[8:8+length_of_header])[0]
        self.number_of_reference_sequences = struct.unpack('<i',first_chunk[8+length_of_header:12+length_of_header])[0]
        rs = length_of_header + 12
        for _ in range(self.number_of_reference_sequences):
            l_name = struct.unpack('<l',first_chunk[rs:rs+4])[0]
            name, l_ref = struct.unpack('<%dsl' % l_name, first_chunk[rs+4:rs+l_name+8]) # disgusting.
            self.chromosome_names.append(name[:-1]) # We dont need the NUL byte.
            self.chromosome_lengths.append(l_ref)
            rs += 8 + l_name
        self.original_binary_header = buffer(first_chunk[:rs])
        self.left_overs = first_chunk[rs:] # bgzip doesnt care about where it drops the division between blocks, so we end up with some
                                           # bytes in the first block that contain read data, not header data.
        for line in self.header_text.split('\n'):
            if line.startswith('@SQ\tSN:'):
                self.chromosomes_from_header.append(line.split('\t')[1][3:]) # God I hope no programs put spaces in the headers instead of tabs..
        if self.chromosomes_from_header != self.chromosome_names:
            print 'ERROR: For some reason, and I have no idea why, the BAM file stores the chromosome name in two locations, the '
            print '       ASCII text header we all know and love, viewable with samtools view -H, and another special binary header'
            print '       which is used to translate the chromosome refID (a number) into a chromosome RNAME when you do bam -> sam.'
            print '       These two headers should always be the same, but apparently someone has messed that up. #YOLO\n'
            print 'Your ASCII header looks like:\n' + self.header_text
            print '\nWhile your binary header has the following chromosomes:'
            print self.chromosome_names
            exit()

        # Now we have the header data stored, but our first call to the generator gave us a block of uncompressed data back that contained way more than
        # just the header data alone. We want this class to be a generator that, on every request, hands back READ data in BAM format uncompressed, starting
        # from the first read, so in order to do that, we need to create a new generator-shim that just regurgitates what was left over from the header parse
        # the first time its called, then uses whatever self.generator would give on all subsequent tries after that.
        # This is how we do that:
        self._iterator = itertools.chain([self.left_overs],self.generator)

    # And this tells python the class is a generator:
    def __iter__(self): return self 
    def next(self): return next(self._iterator)

    def _generator(self):
        try:
            if type(self.file_handle) == str: p = subprocess.Popen(['pigz','-dc',self.file_handle], stdout=subprocess.PIPE)
            elif type(self.file_handle) == file: p = subprocess.Popen(['pigz','-dc'],stdin=self.file_handle, stdout=subprocess.PIPE)
            else: print 'ERROR: I do not know how to open and read from "' + str(self.file_handle) + '"'; exit()
            self.file_handle = p.stdout
            #sys.stderr.write('Using pigz!\n')
        except OSError:
            try:
                if type(self.file_handle) == str:    p = subprocess.Popen(['gzip','--stdout','--decompress','--force',self.file_handle]       , stdout=subprocess.PIPE)
                elif type(self.file_handle) == file: p = subprocess.Popen(['gzip','--stdout','--decompress','--force'], stdin=self.file_handle, stdout=subprocess.PIPE)
                else: print 'ERROR: I do not know how to open and read from "' + str(self.file_handle) + '"'; exit()
                self.file_handle = p.stdout
                #sys.stderr.write('Using gzip!\n')
            except OSError:
                sys.stderr.write('Using internal Python...\n') # We will end up using the python code below. It is faster than the gzip module, but
                                                               # due to how slow Python's zlib module is, it will end up being about 2x slower than pysam.
        data = self.file_handle.read(655360)
        self.bytes_read += 655360
        cache = []
        blocks_left_to_grab = self.blocks_at_a_time
        bs = 0
        checkpoint = 0
        #pool = Pool(processes=3) 
        #def decompress(data): return zlib.decompress(data, 47) # 47 == zlib.MAX_WBITS|32
        decompress = zlib.decompress
        while data:
            if len(data) - bs < 65536:
                data = data[bs:] + self.file_handle.read(35536)
                self.bytes_read += len(data) - bs
                bs = 0

            magic = data[bs:bs+4]
            if not magic: break # a child's heart
            if magic != "\x1f\x8b\x08\x04":
                if magic == 'BAM\1':
                    # The user has passed us already unzipped data, or we're reading from pigz/gzip :)
                    while data:
                        yield data
                        data = self.file_handle.read(35536)
                        self.bytes_read += len(data)
                    raise StopIteration
                elif magic == 'SQLi': print 'OOPS: You have used an SQLite database as your input BAM file!!'; exit()
                else:                 print 'ERROR: The input file is not in a format I understand :('       ; exit()

            try:
                # The gzip format allows compression containers to store metadata about whats inside them. bgzip uses this
                # to encode the virtual file pointers, final decompressed block size, crc checks etc - however we really dont
                # care -- we just want to unzip everything as quickly as possible. So instead of 'following the rules', and parsing this metadata safely,
                # we try to take a short-cut and jump right to the good stuff, and only if that fails we go the long-way-around and parse every bit of metadata properly:

                #cache.append(decompress(data[bs+18:more_bs-8])) ## originally i stored uncompressed data in a list and used map() to decompress in multiple threads. Was not faster.
                #new_zipped_data = data[bs:more_bs]              ## originally i decompressed the data with a standard zlib.decompress(data,47), headers and all. Was slower.
                more_bs = bs + struct.unpack("<H", data[bs+16:bs+18])[0] +1
                cache.append(decompress(data[bs+18:more_bs-8],-15))
                bs = more_bs
            except: ## zlib doesnt have a nice exception for when things go wrong. just "error"
                sys.stderr.write('INFO: Odd bzgip block detected! The author of pybam didnt think this would ever happen... please could you let me know?')
                header_data = magic + data[bs+4:bs+12]
                header_size = 12
                extra_len = struct.unpack("<H", header_data[-2:])[0]
                while header_size-12 < extra_len:
                    header_data += data[bs+12:bs+16]
                    subfield_id = header_data[-4:-2]
                    subfield_len = struct.unpack("<H", header_data[-2:])[0]
                    subfield_data = data[bs+16:bs+16+subfield_len]
                    header_data += subfield_data
                    header_size += subfield_len + 4
                    if subfield_id == 'BC': block_size = struct.unpack("<H", subfield_data)[0]
                raw_data = data[bs+16+subfield_len:bs+16+subfield_len+block_size-extra_len-19]
                crc_data = data[bs+16+subfield_len+block_size-extra_len-19:bs+16+subfield_len+block_size-extra_len-19+8] # I have left the numbers in verbose, because the above try is the optimised code.
                bs = bs+16+subfield_len+block_size-extra_len-19+8
                zipped_data = header_data + raw_data + crc_data
                cache.append(decompress(zipped_data,47)) # 31 works the same as 47.
                # Although the following in the bgzip code from biopython, its not needed if you let zlib decompress the whole zipped_data, header and crc, because it checks anyway (in C land)
                # I've left the manual crc checks in for documentation purposes:
                expected_crc = crc_data[:4]
                expected_size = struct.unpack("<I", crc_data[4:])[0]
                if len(unzipped_data) != expected_size: print 'ERROR: Failed to unpack due to a Type 1 CRC error. Could the BAM be corrupted?'; exit()
                crc = zlib.crc32(unzipped_data)
                if crc < 0: crc = struct.pack("<i", crc)
                else:       crc = struct.pack("<I", crc)
                if expected_crc != crc: print 'ERROR: Failed to unpack due to a Type 2 CRC error. Could the BAM be corrupted?'; exit()

            blocks_left_to_grab -= 1
            if blocks_left_to_grab == 0:
                yield ''.join(cache)
                cache = []
                blocks_left_to_grab = self.blocks_at_a_time
        self.file_handle.close()
        if cache != '': yield ''.join(cache)
