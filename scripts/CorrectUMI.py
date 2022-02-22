#!/usr/bin/env python

import sys
import os
from optparse import OptionParser


###### OPTIONS ######
# options allow for arguments imported from the command line
# they are shown from command line with script.py --help
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to correct UMI reads using UMItools clusterer"

opts.add_option("-b", "--MyBarcodes", help="<Read1> Accepts list of barcodes + UMI in txt")
opts.add_option("-B", "--BClength", default = 16, help="<barcode length")
opts.add_option("-U", "--UMIlength", default = 12, help="<UMI length")
opts.add_option("-o", "--output", help="Output file prefix")

options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()


###### IMPORT MODULES 2 ######
import pandas as pd
from umi_tools import UMIClusterer


if __name__ == "__main__":

	###### INPUT ######
	bf = options.MyBarcodes

	outname = options.output
	bcl = int(options.BClength)
	umil = int(options.UMIlength)
    ##############

	df = pd.read_csv(bf, index_col = None, header = None, names = ['BC_UMI'])
	df['UMI'] = [x[bcl:bcl+umil] for x in df.BC_UMI]         # extract UMIs
	df['BYTEUMI'] = [x.encode() for x in df.UMI]    # Encode
	total = len(df.UMI) - sum(df.UMI=="N"*umil)       # Count valid UMI reads
	UMIcount = df.BYTEUMI.value_counts()            # Count reads for each UMI
	mydict = UMIcount.to_dict()					    

	## Cluster UMIs
	clusterer = UMIClusterer(cluster_method="directional")
	clustered_umis = clusterer(mydict, threshold=1)

	## Dictionary with unique UMIs as key and corrected equivalent as value
	mydict = {}
	for x in clustered_umis:
		correct = x[0]          # correct (more frequent) UMI
		if (len(x) > 1):
			for each in x[1:]:
				mydict[each] = correct
		mydict[correct] = correct
	
	
	df['CORRECTED_UMI'] = [mydict[x] for x in df.BYTEUMI]
	df['CORRECTED_UMI'] = [x.decode() for x in df.CORRECTED_UMI]
	corrected = sum(df.UMI != df.CORRECTED_UMI)
	df.CORRECTED_UMI.to_csv(outname + '.txt', index = False, index_label = None, header = None)

	print("Completed -------------------------------------- \n")
	print("Total UMI reads: " + str(total) + "\n")
	print("Corrected UMI reads: " + str(corrected) + "\n")
	print("% Corrected UMI reads: " + str(round((corrected/total)*100, 2)) + "%\n")
	







