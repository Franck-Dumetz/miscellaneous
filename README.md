# Miscellaneous

Just a few codes that are useful to me and maybe to others

[Picard_SumTable](https://github.com/Franck-Dumetz/miscellaneous/blob/main/Picard_SumTable.py) <br />
Will look for all _ metrics.txt files within a directory and make a CSV file with the percentage of duplication and the estimated size of the library

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
First, put the full path/fasta <br />
second name.csv <br />

[mapping_MultiSample_SlurmArray.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/mapping_MultiSample_SlurmArray.sh) <br />
It is a Slurm array, so all samples will be started at the same time. <br />
Uses Hisat2 to map samples from a generic directory containing many samples, and parses the results with Samtools. <br />

[MapRmDup_noSlurm.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/MapRmDup_noSlurm.sh) <br />
Same as the previous script, but doesn't use Slurm, loops through fastq files from one directory

[Leish_MapRmdup.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/Leish_MapRmdup.sh) <br />
Will map using parameters specific to Leishmania, and other Trypanosomatids (intronic junction deactivated), and remove duplicates using Picard. 



[gfa2fasta.py](https://github.com/Franck-Dumetz/miscellaneous/blob/main/gfa2fasta.py) <br />
Converting a GFA (Graph Fasta) into a Fasta file. <br />
Usage: <br />
    python gfa2fasta.py -i input.gfa -o output.fasta <br />

[RawCounts2CPM.r](https://github.com/Franck-Dumetz/miscellaneous/blob/main/RawCounts2CPM.r) <br />
takes the featureCounts raw count table and transforms it into CPM and logCPM <br />

## submitting ONT basecalled raw signal to SRA
It can only be done so far by converting the POD5 files to FAST5 and basecalling with guppy-4.2.2 (later versions don't do it). <br />
```
pod5 convert to_fast5 /path/to/pod5_dir/*.pod5 --output /path/to/fast5_out/
```
Then use the following script to basecall the fast5 files: [ONT2NCBI.sh](https://github.com/Franck-Dumetz/miscellaneous/blob/main/ONT2NCBI.sh). <br />

[check_fastq_pairing.sh]() (using slurm) will check whether two gzipped FASTQ files are a valid paired-end pair. It confirms the files exist, tests that both gzip files are intact, counts the number of reads in each file, compares the read IDs in order between R1 and R2, and reports the average read length for each mate. All output is written to the text file you provide when you submit the job. <br />
