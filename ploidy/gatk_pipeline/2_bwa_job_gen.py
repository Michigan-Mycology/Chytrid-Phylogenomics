import sys
sys.path.append("/home/amsesk/scripts")
import scriptgen
import numpy as np

NJOBS=5

csv = open(sys.argv[1]).readlines()
csv = np.array(csv)

try:
    jobname_prefix = sys.argv[2]
except:
    jobname_prefix = ""

csv_spl = np.array_split(csv,NJOBS)
for idx,chunk in enumerate(csv_spl):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"{jobname_prefix}_bwa_{idx}",
            account = "tromeara99",
            cpus_per_task = 3,
            mem_per_cpu = 3,
            time = 168
            )
    for line in chunk:
        spl = [x.strip() for x in line.split(',')]
        if len(spl) == 3:
            spl.append("")
        sg.add_command(f"bwa mem -t 3 -a {spl[1]} {spl[2]} {spl[3]} > {spl[0]}.sam")

    sg.write()
