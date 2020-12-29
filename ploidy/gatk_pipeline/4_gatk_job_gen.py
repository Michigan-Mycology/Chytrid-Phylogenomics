import sys
import scriptgen
import numpy as np
import os
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-s", "--samdir", action = "store", required = True, help = "Path to directory containing all of the alignment files.")
parser.add_argument("-a", "--asmdir", action = "store", required = True, help = "Path to the directory containing all of the assemblies.")
parser.add_argument("-d", "--datatable", action = "store", required = True, help = "Path to a 4-column csv that links strain names to the assembly and reads.")
parser.add_argument("--splits", action = "store", required = False, default = 10, help = "Number of distinct batch scripts to split the commands up into. Default: 10")
args = parser.parse_args()

sam_dir = args.samdir
assembly_dir = args.asmdir

lst = [x for x in os.listdir(sam_dir) if x.endswith(".addrg.bam")]

csv_p = [x.split(',') for x in open(args.datatable).readlines()]
csv_p = {x[0]: x[1] for x in csv_p}

lst = np.array_split(lst, int(args.splits) )

for idx,chunk in enumerate(lst):
    sg = scriptgen.SlurmScriptGenerator(
            jobname = f"gatk_{idx}",
            cpus_per_task = 1,
            mem_per_cpu = 25,
            time = 48
            )
    for bam in chunk:
        strain = bam.replace(".dedupped.sorted.addrg.bam", "")
        asm = csv_p[strain]
        asm_path = os.path.join(assembly_dir, asm)

        # Check for required index and dictionary files. If they don't exist, make them.
        faidx_path = f"{asm_path}.fai"
        picard_dict_path = re.sub("[.][a-zA-Z]+$", ".dict", asm_path)

        if not os.path.isfile(faidx_path):
            sg.add_command(f"samtools faidx {asm_path}")

        if not os.path.isfile(picard_dict_path):
            sg.add_command(f"PicardCommandLine CreateSequenceDictionary R={asm_path} O={picard_dict_path}")

        sg.add_command(f"gatk --java-options \"-Xmx24g\" HaplotypeCaller -R {asm_path} -I {bam} -O {bam.replace('.dedupped.sorted.addrg.bam', '.vcf')} -ERC GVCF -A DepthPerAlleleBySample -A MappingQuality -A LikelihoodRankSumTest")

    sg.write()
