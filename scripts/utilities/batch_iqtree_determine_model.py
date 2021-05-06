import sys
sys.path.append(os.path.join(os.environ.get(
    'CHYTRID_PHYLO'), "scripts", "slurm"))
import scriptgen
import numpy as np
import os
import argparse

SUFFIX = "rmgapped"


parser = argparse.ArgumentParser()
parser.add_argument("-a", "--alndir", action = "store", required = True, help = "Path to the directory containing all of the alignments.")
parser.add_argument("-nt", "--numthreads", action="store", required= True, help = "Path to tab-seperated file with alignment paths and numthreads, in that order.")
parser.add_argument("--splits", action = "store", required = False, default = 10, help = "Number of distinct batch scripts to split the commands up into. Default: 10")
args = parser.parse_args()

alndir = args.alndir

numthreads = {}
for line in open(args.numthreads, 'r'):
    spl = [x.strip() for x in line.split('\t')]
    numthreads[spl[0]] = spl[1]

splits = int(args.splits)

file_list = np.array([x for x in os.listdir(alndir) if x.endswith(SUFFIX)])
file_list_spl = np.array_split(file_list, splits)

script_path = os.path.join(os.path.dirname(__file__), "iqtree_determine_model.py")

for idx,chunk in enumerate(file_list_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"iqmf_{idx}",
            cpus_per_task = 12,
            mem_per_cpu = 2,
            time = 48
            )
    for aln in chunk:
        aln_unescape = aln
        aln = scriptgen.escape_pipes_filenames(aln)
        aln_path = os.path.join(os.path.abspath(alndir), aln)
        aln_unescape_path = os.path.join(os.path.abspath(alndir), aln_unescape)
        sg.add_command(f"python {script_path} {aln_path} {numthreads[aln_unescape_path]} > {os.path.basename(aln_path)}.bestmodel")

    sg.write()

