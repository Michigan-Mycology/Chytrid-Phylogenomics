import sys
import os
sys.path.append(os.path.join(os.environ.get(
    'CHYTRID_PHYLO'), "scripts", "slurm"))
import scriptgen
import numpy as np
import re

NJOBS = 10
PEP_PATH = sys.argv[1]
MARKERS_PATH = sys.argv[2]
SUFFIX = "([.]fa$|[.]fasta$)"
SUFFIX_PATTERN = re.compile(SUFFIX)

todo = [os.path.join(PEP_PATH, x) for x in os.listdir(
    PEP_PATH) if re.search(SUFFIX_PATTERN, x) is not None]
todo = np.array(todo)
todo_spl = np.array_split(todo, NJOBS)

jobname = "HmSe"

for idx, chunk in enumerate(todo_spl):
    sg = scriptgen.SlurmScriptGenerator(
        jobname=f"{jobname}_{idx}",
        cpus_per_task=1,
        mem_per_cpu=2,
        time=48
    )
    for line in chunk:
        outline = os.path.basename(line)
        sg.add_command(f"hmmsearch --cpu 1 -E 1e-5 --domtblout {re.sub(SUFFIX_PATTERN, '.odb10.domtbl', outline)} {MARKERS_PATH} {line} > {outline+'.hmmsearch.log'}")

    sg.write()
