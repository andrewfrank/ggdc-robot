#!/usr/bin/env python

# An updated and generalized version of scripts originally developed
# by Katya McGough at American Type Culture Collection

# Controller for automatically submitting jobs to the Genome-to-Genome Distance
# Calculator (GGDC) website: https://ggdc.dsmz.de/ggdc.php, using GGDC v2.1

# todo
# - figure out submission count - why isn't it updating per loop?
# - replace subprocess.call w/ shell = true with safer method

# DEPENDENCIES

import argparse
import os
import sys
import itertools
import numpy
import time
import subprocess

# COMMAND LINE ARGUMENTS

description = ('Automatically submit jobs to the Genome-to-Genome'
               'Distance Calculator (GGDC) website:'
               'https://ggdc.dsmz.de/ggdc.php. Command line arguments must'
               'include either BOTH --queryfile AND --reffile, OR ONLY'
               '--samplefile.')

parser = argparse.ArgumentParser(description = description)
parser.add_argument('--email',
                    help = 'the email address where GGDC will send results')
parser.add_argument('--blastVariant',
                    help = ('the alignment tool used to determine matches '
                    'between query and reference genomes; GGDC recommends '
                    'BLAST+'))
parser.add_argument('--queryfile',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a query '
                    'sequence OR the full path to a fasta file containing '
                    'query sequence(s); identity for the entire list is '
                    'parsed from the first entry'))
parser.add_argument('--reffile',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a reference '
                    'sequence OR the full path to a fasta file containing '
                    'reference sequence(s); identity for the entire list is '
                    'parsed from the first entry'))
parser.add_argument('--samplefile',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for all sample '
                    'sequences OR the full path to a fasta file containing ' 'sample sequences; identity for the entire list is parsed '
                    'from the first entry. WARNING: THIS OPTION IS '
                    'POTENTIALLY EXTREMELY COMPUTATIONALLY INTENSIVE. VERBOSE ' 'FLAG IS RECOMMENDED TO AVOID IMPRACTICAL SUBMISSIONS '))
parser.add_argument('--verbose',
                    help = ('outputs text based checkpoints and interactive '
                    'submission check.'))
args = parser.parse_args()

# FUNCTIONS

# get path this script was called from, useful for creating tmp files
def get_script_path():
    return os.path.dirname(os.path.realpath(sys.argv[0]))

# builds a list of lists for all unique 2 way comparisons from input queryfile
# and reffile
def build_pairs_rq(queryfile, reffile):
    # get line values from query and ref files
    with open(queryfile) as infile:
        qlines = infile.read().splitlines()
    with open(reffile) as infile:
        rlines = infile.read().splitlines()
    # create list of query-reference pairs
    pairs = []
    for qline in qlines:
        for rline in rlines:
            pair = [qline, rline]
            pairs.append(pair)
    # create dictionary where query is a key and value is a list of refs
    pairs_dict = {}
    for pair in pairs:
        if pair[0] in pairs_dict:
            pairs_dict[pair[0]].append(pair[1])
        else:
            pairs_dict[pair[0]] = [pair[1]]
    return(pairs_dict)

# builds a list of lists for all unique 2 way comparisons from input samplefile
# risks making a fuck ton of comparisons, since it's n choose 2
def build_pairs_all(samplefile):
    # get line values from sample file
    with open(samplefile) as infile:
        lines = infile.read().splitlines()
    # create list of query-reference pairs
    pairs = list(itertools.combinations(lines,2))
    # create dictionary where query is a key and value is a list of refs
    pairs_dict = {}
    for pair in pairs:
        if pair[0] in pairs_dict:
            pairs_dict[pair[0]].append(pair[1])
        else:
            pairs_dict[pair[0]] = [pair[1]]
    return(pairs_dict)

# creates qfiles and rfiles for all 2 way comparisons; multiple rfiles of
# roughly the same size per qfile are created when number of refs exceed
# maxrefs value; also outputs useful submission file pair info as a dict
def write_submission_files(pairs_dict, tmp_dir, maxrefs):
    # iterate through query-reference pair dictionary
    files_dict = {}
    for i, (query, refs) in enumerate(pairs_dict.items()):
        # break refs into equal sized chunks 75 lines or smaller
        nchunks = len(refs)/maxrefs + 1
        ref_chunks = numpy.array_split(refs, nchunks)
        # write files
        for j, ref_chunk in enumerate(ref_chunks):
            # write ref file
            rfile_name = os.path.join(tmp_dir,    # create rfile name
                                      'r' + str(i) + '-' +
                                      str(j) + '.txt')
            refs_writeable = '\n'.join(ref_chunk) # format refs
            rfile = open(rfile_name, 'w')         # open rfile
            rfile.write(refs_writeable)           # write to rfile
            # write query file
            qfile_name = os.path.join(tmp_dir,    # create qfile name
                                      'q' + str(i) + '-' +
                                      str(j) + '.txt')
            qfile = open(qfile_name ,'w')         # open qfile
            qfile.write(query)                    # write to qfile
            # write qfile rfile pair info to files_dict
            files_dict[qfile_name] = rfile_name
    return(files_dict)

def submit_ggdc_jobs(crawler_path, files_dict, email, blastVariant):
    job_count = 0
    submission_count = 0
    for qfile, rfile in files_dict.items():
        subprocess.call(['python', crawler_path,
                 email, blastVariant, qfile, rfile],
                shell = True)
        job_count += 1
        submission_count += 1
        print('job count = ' + str(job_count))
        print('submission count = ' + str(submission_count))
        time.sleep(5)
        if submission_count == 6:
            print("6 jobs submitted. Pausing for 25 minutes.")
            time.sleep(1500)
    	    submission_count = 0 # submission count isn't updating???

# SCRIPT

#crawler_path = os.path.join(get_script_path(), 'ggdc-crawler.py')
#tmp_dir = os.path.join(get_script_path(), 'tmp')
crawler_path = 'P:\\Projects\\ggdc-submitter\\src\\ggdc-crawler.py'
tmp_dir = 'P:\\Projects\\ggdc-submitter\\data'

email = 'afrank@atcc.org'
blastVariant = ['GBDP2_BLASTPLUS']
queryfile = 'P:\\Projects\\ggdc-submitter\\data\\query.txt'
reffile = 'P:\\Projects\\ggdc-submitter\\data\\ref.txt'
samplefile = 'P:\\Projects\\ggdc-submitter\\data\\samples.txt'
maxrefs = 75

pairs_dict = build_pairs_rq(queryfile,reffile)
files_dict = write_submission_files(pairs_dict, tmp_dir, maxrefs = 75)
submit_ggdc_jobs(crawler_path, files_dict, email, blastVariant)
