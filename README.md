Automated protein prediction with metaeuk snakemake-ified

Metaeuk is a pretty memory intensive program-- it failed with less than 100 Gb of memory. I think that it might be limited by the read in of the DB? But, not totally sure.

Before using the genemark-es rules you need to activate the environment and run: `cpanm Logger::Simple`
Before using the maker2 environment you need to set up RepeatMasker. Navigage into the conda environment and go to `conda/xxx/share/RepeatMasker` and then run `perl ./configure`. Then follow the instructions. Set the path to the same conda environment but the bin folder.  
