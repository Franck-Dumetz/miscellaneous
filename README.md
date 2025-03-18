# Miscellaneous

Just a few codes that are useful to me and maybe to others

[Summarise Picard output](https://github.com/Franck-Dumetz/miscellaneous/blob/main/Picard_sum.py) <br />

[Summarise many bam files](https://github.com/Franck-Dumetz/miscellaneous/blob/main/summarise_bams.sh) <br />

[ONT read current visualisation](https://github.com/Franck-Dumetz/miscellaneous/blob/main/plot_pod5.py) <br />
code from Kaylee Watson <br />
usage: python plot_pod5.py [pod5_file] [read_ID] <br />

[Batch renaming of directories](https://github.com/Franck-Dumetz/miscellaneous/blob/main/batch_Dir_rename.sh) <br />

[Hisat --new-summary summary table](https://github.com/Franck-Dumetz/miscellaneous/blob/main/Hisat--new-summary_sum.py) <br />
This script creates a table from multiple --new-summary files output by hisat2 <br />

[removing a full sequence in a fasta file using the header](https://github.com/Franck-Dumetz/miscellaneous/blob/main/remove_seq_fasta.py) <br />
Uses the header value without the > sign <br />
usage: <br />
first put the full path/fasta <br />
second name.csv <br />

[mapping_MultiSample_SlurmArray.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/mapping_MultiSample_SlurmArray.sh) <br />
It is slurm array, so all samples will be started at the same time. <br />
Uses Hisat2 to map samples from an generic directory with many smaples and parses the results through Samtools. <br />

[Leish_MapRmdup.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/Leish_MapRmdup.sh) <br />
Will map using parameters specific for Leishmania, and other Trypanosomatids (intronic junction deactivated) and remove duplicates using Picard. 

[gfa2fasta.py](https://github.com/Franck-Dumetz/miscellaneous/blob/main/gfa2fasta.py) <br />
Converting a gfa (graph fasta) into a fasta file. 
Usage:
    python gfa2fasta.py -i input.gfa -o output.fasta
