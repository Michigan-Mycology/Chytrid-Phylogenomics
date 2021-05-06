import sys
import os
sys.path.append(os.path.join(os.environ.get(
    'CHYTRID_PHYLO'), "scripts", "slurm"))
import scriptgen
import numpy as np
import argparse

SUFFIX = "rmgapped"


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alndir", action = "store", required = True, help = "Path to the directory containing all of the alignments.")
parser.add_argument("--splits", action = "store", required = False, default = 10, help = "Number of distinct batch scripts to split the commands up into. Default: 10")
args = parser.parse_args()

alndir = args.alndir
splits = int(args.splits)

file_list = np.array([x for x in os.listdir(alndir) if x.endswith(SUFFIX)])
file_list_spl = np.array_split(file_list, splits)

script_path = os.path.join(os.path.dirname(__file__), "iqtree_determine_numthreads.py")

for idx,chunk in enumerate(file_list_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"iqnt_{idx}",
            cpus_per_task = 12,
            mem_per_cpu = 2,
            time = 48
            )
    for aln in chunk:
        aln = scriptgen.escape_pipes_filenames(aln)
        aln_path = os.path.join(os.path.abspath(alndir), aln)
        sg.add_command(f"python {script_path} {aln_path} > {os.path.basename(aln_path)}.numthreads")

    sg.write()

