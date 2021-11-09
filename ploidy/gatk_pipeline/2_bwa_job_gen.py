import sys
import os
sys.path.append("/home/amsesk/scripts")
import scriptgen
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datatable", action = "store", required = True, help = "Path to a 4-column csv that links strain names to the assembly and reads.")
parser.add_argument("--splits", action = "store", required = False, default = 10, help = "Number of distinct batch scripts to split the commands up into. Default: 10")
args = parser.parse_args()


NJOBS=int(args.splits)

csv = open(args.datatable).readlines()
csv = np.array(csv)

csv_spl = np.array_split(csv,NJOBS)
for idx,chunk in enumerate(csv_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"bwa_{idx}",
            account = "tromeara99",
            cpus_per_task = 3,
            mem_per_cpu = 3,
            time = 168
            )
    for line in chunk:
        spl = [x.strip() for x in line.split(',')]
        if len(spl) == 3:
            spl.append("")

        bwa_index_path = f"{spl[1]}.amb"
        if not os.path.isfile(bwa_index_path):
            sg.add_command(f"bwa index {spl[1]}")

        sg.add_command(f"bwa mem -t 3 -a {spl[1]} {spl[2]} {spl[3]} > {spl[0]}.sam")

    sg.write()
