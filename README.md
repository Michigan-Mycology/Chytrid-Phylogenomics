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
export 
export CHYTRID_PHYLO=/path/to/Chytrid-Phylogenomics/scripts/python
```

3. **Install SCGid** Some of the python scripts included in this repositiory require a FASTA-reading module that is part of SCGid. If you are going to use them, you need to have SCGid installed on your system and be using it in a virtual python3 environment. Doing this will also bring along some other dependencies that are used by these and other scripts. Follow the installation instructions at the [SCGid github repository](https://www.github.com/amsesk/SCGid.git) to install SCGid.

### Before you start...

1. **You need genomes** (specifically predicted proteomes) for the taxa that you want to include in phylogenomic trees. You can get these by sequencing them yourself, from public databases like the [JGI Genome Portal](https://genome.jgi.doe.gov/portal/) or [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/), or most commonly some combination of both. The sequence headers in the FASTAs you get from different places are going to be formatted differently. The first step after downloading all of your predicted proteomes needs to be getting all of those names into the standard NCBI format:
```
LocusTagPrefix|ProteinHeader OptionalDescription
```
There are a variety of ways to do this, but it depends a lot on what the headers in each file are. I don't have any batch scripts for doing this yet, but feel free to ask for help. Some useful command line tools for doing this include `sed` and `awk`.

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

     # Copy the markers over
     # `*.hmm` means all files whose names end with ".hmm"
     cp fungi_odb10/hmms/*.hmm $CHYTRID_PHYLO/markers/.

     # Copy the cutoff scores that BUSCO uses for each marker - we'll use this later
     cp fungi_odb10/scores_cutoff $CHYTRID_PHYLO/odb10_busco_scores_cutoff
     ```

### Pipeline
