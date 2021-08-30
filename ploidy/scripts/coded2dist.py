import sys
import argparse
import pandas as pd
import numpy as np
import logging
from collections import OrderedDict

parser = argparse.ArgumentParser()
parser.add_argument("coded_matrix", action = "store", help = "Path to coded matrix.")
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG)

codes = pd.read_csv(args.coded_matrix, header = 0, index_col = 0)
#print(codes.columns)

#print("Read CSV in.")
#print(codes.shape)

nrow,ncol = codes.shape

distmat = {}

### Remove real shitty samples
row_idx_to_remove = []
for rowname in codes.index:
	row = codes.loc[rowname]
	row_vc = row.value_counts()
	if "N" in row_vc.index:
		ncount = row_vc["N"]
		
		if float(ncount)/float(len(row)) > 0.9:
			row_idx_to_remove.append(rowname)

for rm in row_idx_to_remove:
	logging.warning(f"Removing sample {rm} because it's >90% missing values.")

codes.drop(index=row_idx_to_remove, inplace = True)

nrow,ncol = codes.shape

for p in range(0, nrow):

	logging.info(f"Processing Sample #{p+1}")

	psamp = codes.index[p]
	distmat[psamp] = []
	
	#print(range(p+1, nrow))
	for q in range(0, nrow):

		#print(codes.index[q])
		
		rawdist_fc = 0.0
		comparable_site_count = 0.0
		coliter = 0
		

		for ps,qs in zip(codes.iloc[p,:], codes.iloc[q,:]):
			
			#ps = codes.iloc[p, col]
			#qs = codes.iloc[q, col]

			if ps == "N" or qs == "N":

				# Stop, go to next iteration, and DON'T increase `comparable_site_count`
				continue

			else:
				comparable_site_count += 1.0
				
				if ps == qs:
					pass

				else:
					if ps in ["R", "A"] and qs == "H":
						rawdist_fc += 0.5

					elif ps == "H" and qs in ["R", "A"]:
						rawdist_fc += 0.5

					elif (ps == "A" and qs == "R") or (ps == "R" and qs == "A"):
						rawdist_fc += 1.0

					else:
						raise NotImplementedError("This shouldn't have happened.")

		if (comparable_site_count == 0.0):
			print(psamp, codes.index[q], rawdist_fc, comparable_site_count)

		scaledist_fc = np.true_divide(rawdist_fc, comparable_site_count)

		distmat[psamp].append(str(scaledist_fc))

i = nrow

with open("full_dist.tsv", 'w') as distout:
	distout.write('\t'.join(list([str(x) for x in codes.index])))
	distout.write("\n")
	for key in list(distmat.keys())[::-1]:
		#if i == nrow:
		#	row = ""
		#else:
		
		row = "\t{}".format("\t".join(distmat[key][::-1]))
		
		#zeros = "\t0.0"*i
		distout.write(f"{key}{row}\n")
		i -= 1