#!/usr/bin/env python

###### IMPORT MODULES 1 ######
import os
import re
import regex
import sys
import gzip
import itertools

from optparse import OptionParser

###### OPTIONS ######
# options allow for arguments imported from the command line
# they are shown from command line with script.py --help
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process BC reads and make data suitable for downstream processes"

opts.add_option("-b", "--MyBarcodes", help="<Read1> Accepts list of barcodes in txt")
opts.add_option("-w", "--whitelist", help="<gzip of the universe of valid barcodes to check")
opts.add_option("-l", "--BClength", default = 16, help="<barcode length")
opts.add_option("-d", "--maxdist", help="<Maximun distance from whitelist allowed in order to correct a barcode")

opts.add_option("-n", "--nreads", default = 1000000, help="Number of reads to process in a given chunk")
opts.add_option("-r", "--ncores", default = 4, help="Number of cores for parallel processing.")

opts.add_option("-o", "--output", help="Output file prefix")

options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()


###### IMPORT MODULES 2 ######
from multiprocessing import Pool, freeze_support
from itertools import repeat
from functools import partial

###### INPUT 1 ######
# Read barcodes whitelist
barcodesfilepath = options.whitelist
with gzip.open(barcodesfilepath, "rt") as my_file:
	barcodesR = my_file.readlines()
barcodes = [barcode.rstrip() for barcode in barcodesR]    # remove white spaces
print("Found and imported " + str(len(barcodes)) + " barcodes")	
global barcodes_set 
barcodes_set = set(barcodes)

# Dictionary of possible substitutions (ALPHABET_MINUS)
DNA_ALPHABET = 'AGCT'
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}
ALPHABET_MINUS['N'] = set(DNA_ALPHABET)   
MAXDIST_CORRECT = 2


###### FUNCTIONS ######
def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


def gen_nearby_seqs(seq,wl_idxs=None, maxdist=1):
	"""
	Generate all sequences with at most maxdist changes from seq that are in a provided whitelist
	Returns only the modified sequences that are in whitelist
	Output is a generator. It can be printed with a loop
	If the input sequence is in whitelist => empty output
	If distance with any sequence in barcode is > maxdist => empty output
	** wl_idxs is not used
    """	
	allowed_indices = [i for i in range(len(seq)) if seq[i] != 'N']      ## count no Ns
	required_indices = tuple([i for i in range(len(seq)) if seq[i] == 'N'])     ## count Ns
	mindist = len(required_indices)
	if mindist > maxdist:     
		return
	#if mindist == 0:
	#	mindist = mindist + 1
	for dist in range(mindist, maxdist + 1):    # Iterate over the allowed distances (+1)
		for modified_indices in itertools.combinations(allowed_indices, dist - mindist):   # Iterate over possible combinations of substitutions
			indices = set(modified_indices + required_indices)                             # Add index of required substitution, when there's an N

			for substitutions in itertools.product(*[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
				new_seq = ''.join(substitutions)
				if new_seq in barcodes_set:        ## barcodes_set = whitelist
					yield new_seq
					
#------ 

def correct_barcode(seq, maxdist=1, mylength=16):
	"""
	Correct barcodes not in whitelist
    """
	if seq in barcodes_set:
		return seq

	for test_str in gen_nearby_seqs(seq, maxdist=maxdist):
		return(test_str)

	return "N"*mylength


###########################################################################

if __name__ == "__main__":

	###### INPUT 2 ######
	bf = options.MyBarcodes
	bcl = int(options.BClength)
	outname = options.output
	outfilename = outname + ".txt"
	cpu = int(options.ncores)
	n = int(options.nreads)
	maxdist = int(options.maxdist)

## how to read BC file, iterate over it and write the corrected version
## how to make it in parallel
if __name__ == "__main__":
	discarded = 0                 ## BC correction counters
	changed = 0
	total = 0
	with open(bf, "rt") as f1:
		# Establish iterators
		it1 = batch_iterator(f1, n)

		# iterate over batches of length n barcodes
		with open(outfilename, 'w') as out_write:
			for i, batch1 in enumerate(it1):
				batch1 = [x.strip() for x in batch1]   # remove \n character
				batch1 = [x[0:bcl] for x in batch1]   # select barcode
				nbcs = len(batch1)
				print("Reading and correcting " + str(nbcs) + " barcodes")
				total+=nbcs
				# parallel process the barcode processing and accounting of failures.
				pool = Pool(processes=cpu)
				pm = pool.map(partial(correct_barcode, maxdist = maxdist, mylength = bcl), batch1)
				pool.close()
				for bc, bc_corrected in zip(batch1, pm):     # keep count of discarded and corrected barcodes
					if bc != bc_corrected:
						changed+=1
						if bc_corrected == "N"*bcl:
							discarded+=1
				pm = [x + "\n" for x in pm]           # add \n to write into file
				out_write.writelines(pm)

	out_write.close()
	f1.close()

	corrected = changed - discarded

	print("Completed -------------------------------------- \n")
	print("Total barcode reads: " + str(total) + "\n")
	print("Corrected barcode reads: " + str(corrected) + "\n")
	print("% Corrected barcode reads: " + str(round((corrected/total)*100, 2)) + "%\n")
	print("Discarded barcode reads: " + str(discarded) + "\n")
	print("% Discarded barcode reads: " + str(round((discarded/total)*100, 2)) + "%\n")






















