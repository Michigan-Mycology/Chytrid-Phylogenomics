# Chytrid-Phylogenomics
Repository of scripts and shared data for chytrid phylogenomics analyses

### Set some environmental variables

If you're going to be using the automated batch gene tree filtering algorithm, you need to set an environmental variable to tell python where some python library files are located. You can set this in your current session by simply entering the following command:

```
CHYTRID_PHYLO_PY=/path/to/Chytrid-Phylogenomics/scripts/python
```

You can have this set automatically everytime you login by adding the following line to `/home/uniqname/.bashrc`:

``` /home/uniqname/.bashrc
export CHYTRID_PHYLO_PY=/path/to/Chytrid-Phylogenomics/scripts/python
```
### Pipeline
