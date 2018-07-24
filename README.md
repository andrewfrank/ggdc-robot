# ggdc-robot v0.0.4

Automatically submit jobs to the Genome-to-Genome Distance Calculator (GGDC)
website: https://ggdc.dsmz.de/ggdc.php, using GGDC v2.1.

## Prerequisites

ggdc-robot relies on the following Python packages:
- itertools
- numpy
- requests

## Installing

Move ggdc-robot.py to a directory of your choice. Run using the command:
```
python ggdc-robot.py
```

## Running

When running, ggdc-robot.py creates the directory "submissions" in the same directory where the script is run from. This directory stores files to submit to GGDC.

## Usage
```
usage: ggdc-robot.py [-h] (--samplefile FILE | --queryfile FILE)
                     [--reffile FILE] --email EMAIL ADDRESS
                     (--slotusage PERCENT | --bruteforce)
                     [--blastVariant {GBDP2_BLASTPLUS,GBDP2_BLAT,GBDP2_BLASTZ,GBDP2_WU-BLAST,GBDP2_MUMMER}]
                     [--timedwait MINUTES JOBS] [--version]

arguments:
 -h, --help            show this help message and exit
 --samplefile FILE, -s FILE
                       the full path to a text file where each new line
                       contains EITHER the NCBI accession number for all
                       sample sequences OR the full path to a fasta file
                       containing sample sequences; identity for the entire
                       list is parsed from the first entry; NOTICE: THIS
                       OPTION IS POTENTIALLY EXTREMELY COMPUTATIONALLY
                       INTENSIVE.
 --queryfile FILE, -q FILE
                       the full path to a text file where each new line
                       contains EITHER the NCBI accession number for a query
                       sequence OR the full path to a fasta file containing
                       query sequence(s); identity for the entire list is
                       parsed from the first entry.
 --reffile FILE, -r FILE
                       the full path to a text file where each new line
                       contains EITHER the NCBI accession number for a
                       reference sequence OR the full path to a fasta file
                       containing reference sequence(s); identity for the
                       entire list is parsed from the first entry.
 --email EMAIL ADDRESS, -e EMAIL ADDRESS
                       the email address where GGDC will send results.
 --slotusage PERCENT, -u PERCENT
                       enable slot usage waiting mode; this mode forces ggdc-
                       robot to pause for 10 minutes before attempting to
                       submit another job when GGDC servers reach this
                       specificied capacity. E.g. --slotusage 50 will prompt
                       ggdc-robot to wait 10 minutes when GGDC server slot
                       usage reaches 50 percent.
 --bruteforce, -f      enable brute force mode; this mode forces ggdc-robot
                       to submit jobs to the GGDC server even when the server
                       load is at 100 percent; ATTENTION: jobs may fail due
                       to GGDC job queue limits.
 --blastVariant {GBDP2_BLASTPLUS,GBDP2_BLAT,GBDP2_BLASTZ,GBDP2_WU-BLAST,GBDP2_MUMMER},
 -b {GBDP2_BLASTPLUS,GBDP2_BLAT,GBDP2_BLASTZ,GBDP2_WU-BLAST,GBDP2_MUMMER}
                       the alignment tool used to determine matches between
                       query and reference genomes; GGDC recommends
                       GBDP2_BLASTPLUS.
 --timedwait MINUTES JOBS, -t MINUTES JOBS
                       enable timed waiting mode; this mode forces ggdc-robot
                       to wait the X minutes every Y jobs submitted. E.g.
                       --wait 25 6 will force ggdc-robot to wait 25 minutes
                       between every set of 6 jobs submitted. Can be combined
                       with --usagewait if desired.
 --version             show program's version number and exit
```
