import sys
import os

L50_HEADER_DIR = "/home/aimzez/DATA/pursuit/ploidy_final/assembly_l50_headers"
DEPTH_DIR = "/home/aimzez/DATA/pursuit/ploidy_final/contig.depth"

L50_HEADER_SUFF = "l50_headers"
DEPTH_SUFF = "sorted.bam.depth.ContigMean"

prefix_to_asm = {}
with open(sys.argv[1], 'r') as datatable:
    for line in datatable:
        spl = line.split(',')
        prefix_to_asm[spl[0]] = os.path.basename(spl[1])

for depth_file in [x for x in os.listdir(DEPTH_DIR) if x.endswith(DEPTH_SUFF)]:
    prefix = depth_file.replace(f".{DEPTH_SUFF}", "")
    
    l50_header_path = os.path.join(L50_HEADER_DIR, f"{prefix_to_asm[prefix]}.{L50_HEADER_SUFF}")
    l50_headers = [x.strip() for x in open(l50_header_path, 'r').readlines()]
    l50_headers = [x.split(" ")[0] for x in l50_headers]

    depths_path = os.path.join(DEPTH_DIR, depth_file)
    depths = [x.strip() for x in open(depths_path, 'r').readlines()]

    first_line = depths.pop(0)

    with open(f"{prefix}.{DEPTH_SUFF}.L50", 'w') as out:
        out.write(first_line)
        out.write("\n")
        for line in depths:
            spl = line.split("\t")
            if spl[0] in l50_headers:
                out.write(f"{line.strip()}\n")
            


