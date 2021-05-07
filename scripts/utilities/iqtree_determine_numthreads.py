import sys
import subprocess
import re
import time
import os

pat = re.compile("BEST[ ]NUMBER[ ]OF[ ]THREADS[:][ ]([0-9]+)")

alignment = sys.argv[1]
logfile = f"{os.path.basename(alignment)}.log"
cmd = [
        "iqtree",
        "-m", "MF",
        "-s", f"{alignment}",
        "-T", "AUTO",
        "--prefix", f"{os.path.basename(alignment)}"
        ]

p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

out,err = p.communicate()
print(err)

sys.exit()
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




