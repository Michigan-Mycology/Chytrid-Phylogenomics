import sys
sys.path.append("/home/amsesk/scripts/")
import scriptgen
import numpy as np
import os

NJOBS=2
UNALN_PATH = os.path.abspath(sys.argv[1])
MARKERS_PATH = os.path.abspath(sys.argv[2])

todo = [x for x in os.listdir(UNALN_PATH) if x.endswith(".fasta")]
todo = np.array(todo)
todo_spl = np.array_split(todo, NJOBS)

jobname = "HmAl"

for idx,chunk in enumerate(todo_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"{jobname}_{idx}",
            cpus_per_task = 1,
            mem_per_cpu = 2,
            time = 48
            )
    for line in chunk:
        sg.add_command(f"hmmalign --trim --amino -o {line.replace('.fasta','.aa.msa')} {os.path.join(MARKERS_PATH, line.replace('.fasta','.hmm'))} {line}")
        sg.add_command(f"esl-reformat --replace=x:- --gapsym=- -o {line.replace('.fasta','.esltmp')} afa {line.replace('.fasta','.aa.msa')}")
        sg.add_command("perl -p -e 'if (! /^>/) {{ s/[ZBzbXx\*]/-/g }}' {} > {}".format(line.replace('.fasta','.esltmp'), line.replace('.fasta','.clnaln')))
        sg.add_command(f"rm {line.replace('.fasta','.esltmp')}")
        sg.add_command(f"trimal -resoverlap 0.50 -seqoverlap 60 -in {line.replace('.fasta','.clnaln')} -out {line.replace('.fasta','.aa.filter')}")
        sg.add_command(f"trimal -fasta -in {line.replace('.fasta','.aa.filter')} -out {line.replace('.fasta','.aa.trim')}")

    sg.write()
