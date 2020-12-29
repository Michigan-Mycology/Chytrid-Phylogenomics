import sys
sys.path.append(os.path.join(os.environ.get('CHYTRID_PHYLO'), "scripts", "slurm"))
import scriptgen
import numpy as np
import os

NJOBS=10
PEP_PATH = sys.argv[1]
MARKERS_PATH = sys.argv[2]

todo = [os.path.join(PEP_PATH, x) for x in os.listdir(PEP_PATH) if x.endswith(".aa.fasta")]
todo = np.array(todo)
todo_spl = np.array_split(todo, NJOBS)

jobname = "HmSe"

for idx,chunk in enumerate(todo_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"{jobname}_{idx}",
            cpus_per_task = 1,
            mem_per_cpu = 2,
            time = 48
            )
    for line in chunk:
        sg.add_command(f"hmmsearch --cpu 1 -E 1e-5 --domtblout {line.replace('.aa.fasta', 'odb10.domtbl')} {MARKERS_PATH} {line} > {line+'.hmmsearch.log'}")

    sg.write()
