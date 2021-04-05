import sys
sys.path.append(os.path.join(os.environ.get(
    'CHYTRID_PHYLO'), "scripts", "slurm"))
import scriptgen
import numpy as np
import os

NJOBS = 20
ALN_PATH = os.path.abspath(sys.argv[1])

todo = [f"{os.path.join(ALN_PATH,x)}" for x in os.listdir(ALN_PATH) if x.endswith(".aa.trim.rmgapped")]
todo = np.array(todo)
todo_spl = np.array_split(todo, NJOBS)

jobname = "FaTr"

for idx, chunk in enumerate(todo_spl):
    sg = scriptgen.SlurmScriptGenerator(
        jobname=f"{jobname}_{idx}",
        cpus_per_task=1,
        mem_per_cpu=2,
        time=1
    )
    sg.add_command("module load singularity")
    sg.add_command("module load fasttree")
    for line in chunk:
        outbase = os.path.basename(line)
        sg.add_command(f"fasttree -gamma < {line} > {outbase.replace('.trim.rmgapped','.tre')}")

    sg.write()
