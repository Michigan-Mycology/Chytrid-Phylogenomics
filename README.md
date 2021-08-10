# Chytrid-Phylogenomics
Repository of scripts and shared data for chytrid phylogenomics analyses. 

An emphasis of the workflow described here is that arbitrarily choosing the "best hits" for a phlogenomic marker may not always produce the best gene trees, negatively impacting downstream concatenated analyses and so on.

### Setup
1. **Clone this repository**
```
git clone https://github.com/Michigan-Mycology/Chytrid-Phylogenomics.git
```

2. **Set an environmental variable that points to where you have this repository downloaded.** Do this by adding the following line to the the file at `/home/uniqname/.bashrc` or just typing it into your current session. Adding it to `/home/uniqname/.bashrc` will make this permanent.

``` 
export CHYTRID_PHYLO=/path/to/Chytrid-Phylogenomics/scripts/python
```

3. **Install SCGid** Some of the python scripts included in this repositiory require a FASTA-reading module that is part of SCGid. If you are going to use them, you need to have SCGid installed on your system and be using it in a virtual python3 environment. Doing this will also bring along some other dependencies that are used by these and other scripts. 
     
     Start to follow the installation instructions at the [SCGid github repository](https://www.github.com/amsesk/SCGid.git) to install SCGid, but you can stop after you enter:
     
     `python setup.py develop`

### Before you start...

1. **You need genomes** (specifically predicted proteomes) for the taxa that you want to include in phylogenomic trees. You can get these by sequencing them yourself, from public databases like the [JGI Genome Portal](https://genome.jgi.doe.gov/portal/) or [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/), or most commonly some combination of both. The sequence headers in the FASTAs you get from different places are going to be formatted differently. The first step after downloading all of your predicted proteomes needs to be getting all of those names into the standard NCBI format:
     ```
     LocusTagPrefix|ProteinHeader OptionalDescription
     ```
     There are a variety of ways to do this, but it depends a lot on what the headers in each file are. I don't have any batch scripts for doing this yet, but feel free to ask for help. Some useful command line tools for doing this include `sed` and `awk`. Here are some examples working with proteomes from [GenBank](https://github.com/Michigan-Mycology/Lab-Code-and-Hacks/tree/master/Phylogenomics/processing_genbank_files) and [JGI](https://github.com/Michigan-Mycology/Lab-Code-and-Hacks/tree/master/Phylogenomics/processing_jgi_files).
     
     Put the final predicted proteome FASTA files at `$CHYTRID_PHYLO/data`. You will have to make this directory first, `mkdir $CHYTRID_PHYLO/data`.

2. **You need phylogenomic markers** that you want to pull out of the predicted proteomes for taxa you want to include in phylogenomic trees. For this particular project we opted to use the markers used by [BUSCO](https://busco.ezlab.org/) and contained in the `fungi_odb10` database. There are other varieties of BUSCO markers that are tuned for different groups of organisms.
     - You can find fungi_odb10 [HERE](https://busco-data.ezlab.org/v4/data/lineages/fungi_odb10.2019-12-13.tar.gz)
     - You can find a list of all the odb10 databases [HERE](https://busco.ezlab.org/busco_v4_data.html)

     You can download `fungi_odb10` to your computer from the terminal by entering
     ```
     # Download from the internet
     wget https://busco-data.ezlab.org/v4/data/lineages/fungi_odb10.2019-12-13.tar.gz

     # Untar/Unzip compressed archive
     tar xvf fungi_odb10.2019-12-13.tar.gz

     # Now you should have a directory called `fungi_odb10`
     ```

     The odb10 archive comes with several directories and files. We don't need all of them to run this pipeline, so I like to copy the important stuff (for us) out of that archive and put them in a working directory that I'll use for the entire pipeline.
     ```
     # Make a direcotry for the markers
     mkdir $CHYTRID_PHYLO/markers
     mkdir $CHYTRID_PHYLO/markers/hmm

     # Copy the markers over
     # `*.hmm` means all files whose names end with ".hmm"
     cp fungi_odb10/hmms/*.hmm $CHYTRID_PHYLO/markers/hmm/.

     # Copy the cutoff scores that BUSCO uses for each marker - we'll use this later
     cp fungi_odb10/scores_cutoff $CHYTRID_PHYLO/odb10_busco_scores_cutoff
     ```

### Pipeline

Make sure you have markers and genomes downloaded and FASTA headers reformatted as discussed above before continuing.

1. Combine all the hmm markers into a single file and format as a database.

     Install/load hmmer. If you are on greatlakes, hmmer is already installed, but you need to load it:
     ```
     module load Bioinformatics hmmer
     ```
     Combine files and build database like this:
     ```
     cd $CHYTRID_PHYLO/markers/
     cat hmm/*.hmm > odb10_combined.hmm
     hmmpress odb10_combined.hmm
     ```
     Now you should have `.h3p`, `.h3m`, `.h3f`, and `h3i` files in addition to `odb10_combined.hmm`.
     
2. Search predicted proteomes with marker hmm models to get hits in domtbl format.
     
     The essential command that needs to be called for each proteome is:
     ```
     hmmsearch --cpu 1 -E 1e-5 --domtblout predicted_proteome.domtbl  /path/to/odb10_combined.hmm /path/to/predicted_proteome.fasta
     ```
     
     **If you are working on greatlakes through the James Lab** there is a python helper script you can use to autogenerate slurm scripts to run this command for all your FASTA files as long as they're stored together in a single directory.
     ```
     python $CHYTRID_PHYLO/scripts/slurm/1_hmmsearch.py /path/to/FASTA/directory /path/to/odb10_combined.hmm
     ```
     
     Now you should have 10 HmSe_X.sh files which can each be submitted to slurm, the results of which should be collected in a separate directory:
     ```
     mkdir $CHYTRID_PHYLO/search
     cd $CHYTRID_PHYLO/search
     for i in $(ls /path/to/slurm/scripts/HmSe_*.sh); sbatch $i
     ```
     

     
