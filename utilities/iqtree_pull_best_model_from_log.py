import sys
import os
import re
import argparse
import pandas as pd

pattern = re.compile("^Best[-]fit[ ]model[:][ ]([^ ]+).*$")
parser = argparse.ArgumentParser()
parser.add_argument("--logdir", action="store", help="Path to directory containing IQTree ModelFinder logfiles.")
args = parser.parse_args()

ldict = []
for logfile in [x for x in os.listdir(args.logdir) if x.endswith(".log")]:
    file_base = logfile.split('.')[0]
    for line in open(os.path.join(args.logdir, logfile)):
        s = re.match(pattern, line)
        if s is not None:
            ldict.append({'marker': file_base, 'best_model': s.group(1)})

df = pd.DataFrame(ldict)
df.to_csv("iqtree_marker_by_best_models.tsv", sep="\t", index=False)

popularity = df.best_model.value_counts().reset_index()
popularity.columns = ["model", "n_markers_best_for"]
popularity.to_csv("iqtree_model_popularity.tsv", sep="\t", index=False)
