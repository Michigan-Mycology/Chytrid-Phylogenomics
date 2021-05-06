import sys
import subprocess
import re
import time
import os

pat = re.compile("Best-fit[ ]model[:][ ]([a-zA-Z0-9+]+)[ ].*")

alignment = sys.argv[1]
nthreads = sys.argv[2]

logfile = f"{os.path.basename(alignment)}.log"
cmd = [
        "iqtree",
        "-m", "MF",
        "-s", f"{alignment}",
        "-T", f"{nthreads}",
        "--prefix", f"{os.path.basename(alignment)}"
        ]

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

while True:
    time.sleep(5)
    with open(logfile, 'r') as lf:
        for line in lf.readlines():
            s = re.search(pat, line)
            if s is not None:
                print (f"{alignment}\t{s.group(1)}")
                while p.poll() is None:
                    p.terminate()

                break

    if p.poll() is not None:
        break



