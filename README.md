# LinkGe

[![DOI](https://zenodo.org/badge/123921192.svg)](https://zenodo.org/badge/latestdoi/123921192)

LinkGe is a perl program designed to quickly generate linkage polymorphism frequency reports from a sequence alignment given in bam format. It works in both nucleotide and amino acid space. Additionally it works with paired and non-paired data, treating paired reads as one read when searching for linkages. No theoretical limit on the amount of the positions to be queried. It takes under 1 second to process a position with over 3000 reads of coverage.

LinkGe is a perl script designed to take in an Illumina alignment in BAM format and output linkage polymorphism (LP) frequencies. LPs are determined by user-defined nucleotide positions from the command line or in the PARAM file. At each of those positions LinkGe collects the nucleotide base call for each read and then collects that data

Options:
```
-v    --version     print current version of LinkGe
-r    --range       allow an input range
-o    --out         define the prefix for the output file
-b    --batch       process all of the files in a given folder
-h    --help        show this screen
-t    --threaded    run each file from a batch folder on a separate thread 
-F    --FASTA       export FASTA files for each linkage polymorphism
-f    --file        define positions from a csv or tsv file
-l    --list        define positions from a command line list
-A    --Analysis    amino acid analysis
```
Usage:
```
$ LinkGe.pl [OPTIONS] [PARAM file] [alignment.bam]
```
Dependencies:
```
Bio::DB::Sam
BioPerl
Samtools version 0.1.15
```
Aquire Bio::DB::Sam from CPAN and install all associated dependencies.

Example of using a range to define positions
```
$ perl LinkGe.pl -r 2189-2275 alignment.bam
```
Example of using a list to define positions
```
$ perl LinkGe.pl  -l 1190,1201,1205,1224,1232 alignment.bam
```
Examples of accepted PARAM file format:
```
$ perl LinkGe.pl  -f /path/to/your/file.txt alignment.bam
```
Comma separated:
```
position_1,8
position_2,25
position_3,352
position_4,700
...
position_n,N
```
OR

Tab separated:
```
position_1  8
position_2  25
position_3  352
position_4  700
...
position_n  N
```
Both formats must have Unix end-line characters (LF).

LinkGe was first published in:
P.R. Wilker, J.M. Dinis, G. Starrett, M. Imai, M. Hatta, C.W. Nelson, D.H. O’Connor, et al. “Selection on Haemagglutinin Imposes a Bottleneck during Mammalian Transmission of Reassortant H5N1 Influenza Viruses.” Nature Communications 4 (2013). https://doi.org/10.1038/ncomms3636. 

Documentation
```
0.1.0
-First "working" version

0.1.1
-allows for command line defined location for bam file
-prints usage screen for invalid number of arguments

0.1.2
-Writes one file, linkage_frequencies.txt to whichever directory the user is in when starting PEGen.

0.1.3
-Using fetch instead of pileup to collect read name and base calling information (pileup was missing up to 12,000 reads at certain locations. That's why there was the undef problem).
-Fetch also makes it magnitudes faster.
-Deletions are now reported as a "-" in the final linkage report.
-Put new comparison commands into a subroutine so now the whole program is scaleable. Minimum two positions must be given, but I have not tested going beyond 12 positions.
-Can now accept the PARAM file in both tab separated and comma separated value formats.
-Calculates frequency in output file.
-Changed pattern reporting format from e.g. ATGTAG to A107:T122:G130:T138:A164:G168:

0.1.4
-first implementation of options: version, batch, out, and help
-default output filename format is now bamfilename_timestamp_year_month_day.txt
-suppressed uninitialized and substring warnings to prevent warnings when occasionally it looks outside the length of a sequence string
-aesthetic changes to usage string
-separated version, author, etc into separate scalar variable to more easily update.
-small rearrangements to code to allow option functionality to work
-realized that fasta file is not required
-modifications to code to allow for non paired-end data

0.2.0
-option for multithreaded batch processing
-option for exporting FASTA files
-a few error checking steps for conflicting options
-numbering system to associate FASTA files with linkage polymorphisms in report text file
-renamed program to LinkGe

0.2.1
-changed the position input options. No longer is the user constrained to use the external param file to define position information
-Addition of the Analysis option (explained elsewhere)

0.2.2
-small bug fixes and moved onto github
```
