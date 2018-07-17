#!/usr/bin/env python

__version__ = '0.0.3'

# Controller for automatically submitting jobs to the Genome-to-Genome Distance
# Calculator (GGDC) website: https://ggdc.dsmz.de/ggdc.php, using GGDC v2.1.
# Updated and generalized version of scripts originally developed by Katya
# McGough at American Type Culture Collection.

# DEPENDENCIES

import argparse
import os
import sys
import itertools
import numpy
import time
import mechanize
from bs4 import BeautifulSoup

# COMMAND LINE ARGUMENTS

description = ('Automatically submit jobs to the Genome-to-Genome'
               'Distance Calculator (GGDC) website:'
               'https://ggdc.dsmz.de/ggdc.php. Command line arguments must '
               'include either BOTH --queryfile AND --reffile, OR ONLY'
               '--samplefile.')

parser = argparse.ArgumentParser(description = description)
parser.add_argument('--email','-e',
                    help = 'the email address where GGDC will send results',
                    required = True)
parser.add_argument('--blastVariant','-b',
                    help = ('the alignment tool used to determine matches '
                    'between query and reference genomes; GGDC recommends '
                    'GBDP2_BLASTPLUS'),
                    choices = ['GBDP2_BLASTPLUS','GBDP2_BLAT','GBDP2_BLASTZ',
                    'GBDP2_WU-BLAST','GBDP2_MUMMER'],
                    default = 'GBDP2_BLASTPLUS')
ifiles = parser.add_mutually_exclusive_group(required = True)
ifiles.add_argument('--samplefile','-s',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for all sample '
                    'sequences OR the full path to a fasta file containing '
                    'sample sequences; identity for the entire list is parsed '
                    'from the first entry; WARNING: THIS OPTION IS '
                    'POTENTIALLY EXTREMELY COMPUTATIONALLY INTENSIVE. VERBOSE '
                    'FLAG IS RECOMMENDED TO AVOID IMPRACTICAL SUBMISSIONS '))
ifiles.add_argument('--queryfile','-q',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a query '
                    'sequence OR the full path to a fasta file containing '
                    'query sequence(s); identity for the entire list is '
                    'parsed from the first entry'))
parser.add_argument('--reffile','-r',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a reference '
                    'sequence OR the full path to a fasta file containing '
                    'reference sequence(s); identity for the entire list is '
                    'parsed from the first entry'),
                    required = '--queryfile' in sys.argv)
parser.add_argument('--bruteforce', '-f',
                    help = ('enable brute force mode; this mode forces '
                    'ggdc-robot to submit jobs to the GGDC server even when the'
                    'server load is at 100%; ATTENTION: jobs may fail due to '
                    'GGDC job queue limits'),
                    action = 'store_true')
parser.add_argument('--wait', '-w',
                    help = ('enable waiting mode; this mode forces '
                    'ggdc-robot to wait the X minutes every Y jobs '
                    'submitted. E.g. --wait 25 6 will force ggdc-robot to '
                    'wait 25 minutes between every set of 6 jobs submitted.'),
                    type = int,
                    nargs = 2)
parser.add_argument('--verbose','-v',
                    help = ('outputs text based checkpoints and interactive '
                    'submission check.'),
                    action = 'store_true')
parser.add_argument('--version',
                    action = 'version',
                    version = __version__)
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
def write_submission_files(pairs_dict, submissions_dir, maxrefs):
    # iterate through query-reference pair dictionary
    files_dict = {}
    for i, (query, refs) in enumerate(pairs_dict.items()):
        # break refs into equal sized chunks 75 lines or smaller
        nchunks = len(refs)/maxrefs + 1
        ref_chunks = numpy.array_split(refs, nchunks)
        # write files
        for j, ref_chunk in enumerate(ref_chunks):
            # write ref file
            rfile_name = os.path.join(submissions_dir,    # create rfile name
                                      'r' + str(i) + '-' +
                                      str(j) + '.txt')
            refs_writeable = '\n'.join(ref_chunk) # format refs
            rfile = open(rfile_name, 'w')         # open rfile
            rfile.write(refs_writeable)           # write to rfile
            # write query file
            qfile_name = os.path.join(submissions_dir,    # create qfile name
                                      'q' + str(i) + '-' +
                                      str(j) + '.txt')
            qfile = open(qfile_name ,'w')         # open qfile
            qfile.write(query)                    # write to qfile
            # write qfile rfile pair info to files_dict
            files_dict[qfile_name] = rfile_name
    return(files_dict)

def check_submission_format(file):
    if os.path.isfile(file):
        with open(file) as infile:
            first_line = infile.readline().rstrip()
        if os.path.isfile(first_line):
            return('filepath')
        else:
            return('accession')
    else:
        sys.exit(file + ' not found.')

def get_ggdc_status(url):
    # browse to GGDC website
    br = mechanize.Browser()
    page = br.open(url)

    # check server load
    page_html = BeautifulSoup(page, 'html.parser')
    status_html = page_html.find('div', {'class': 'progress'})
    status = status_html.get_text()
    status_message = 'Current GGDC server load:' + status
    return(status_message)

def submit_file_values(file, type, format):
    with open(file) as infile:
        lines = infile.read().splitlines()

    if format == 'accession':
        if type == 'query':
            control = 'targetName'
            sep = ' '
        elif type == 'ref':
            control = 'refGenbank'
            sep = '\r\n'
        submission = sep.join(lines)
        try:
            br.form.set_value(submission, control)
            #print(submission)
        except:
            print('Error submitting ' + submission + ' ' + type + '.')

    elif format == 'filepath':
        if type == 'query': name = 'targetGenome'
        elif type == 'ref': name = 'multipleRefGenomes[]'
        for submission in lines:
            # ----- START ERROR FOR MULTIPLE REF FASTAS ----- #
            if len(lines) > 1:
                sys.exit('Sorry, submission of multiple reference fasta files '
                         'is not yet supported. Exiting.')
            else:
            # ----- END ERROR FOR MULTIPLE REF FASTAS ----- #
                try:
                    br.form.add_file(open(submission),
                                     'text/plain',
                                     submission,
                                     name = name)
                except:
                    print('Error submitting ' + submission + ' ' + type + '.')

    else: sys.exit('Unable to submit ' + type + ' ' + format + '. Exiting.')

def ggdc_submit(url, email, blastVariant, queryfile, reffile):
    # browse to GGDC website
    br = mechanize.Browser()
    br.open(url)

    # begin filling out GGDC form
    br.select_form('Form')                          # GGDC form name is "Form"
    br.form.set_value(email, 'email')               # fill in email form
    br.form.set_value([blastVariant], 'blastVariant') # fill in BLAST method form

    # fill in query form per format submitted
    queryformat = check_submission_format(queryfile)
    submit_file_values(file = queryfile,
                       type = 'query',
                       format = queryformat)

    # fill in ref form per format submitted
    refformat = check_submission_format(reffile)
    submit_file_values(file = reffile,
                       type = 'ref',
                       format = refformat)

    # submit GGDC job
    submission = br.submit()

    # get GGDC job response
    submission_html = BeautifulSoup(submission, 'html.parser')
    try:
        response_html = submission_html.find('div', {'class': 'alert alert-success'})
        response = response_html.get_text()
    except:
        response_html = submission_html.find('div', {'class': 'panel-body alert-danger'})
        response = response_html.get_text()
    print(response)

# iteratively submits each qfile-rfile pair to GGDC using ggdc-crawler.py;
# currently pauses for 25 minutes every 6th submission
def submit_ggdc_jobs(url, email, blastVariant, files_dict, bruteforce, wait):

    submission_count = 0
    jobs_requested = len(files_dict)
    print('Jobs requested = ' + str(jobs_requested))

    for job_count, (qfile, rfile) in enumerate(files_dict.items()):
        status = get_ggdc_status(url)
        print(status)

        if bruteforce is not None:
           while '100%' in status:
               print(('All GGDC server slots are currently used. Waiting '
                      '10 minutes before attempting additional submissions.'))
               print('Waiting to submit job ' + str(job_count) + '.')
               if wait is not None:
                   print('This is job ' + str(submission_count) +
                         ' of this submission set.')
               time.sleep(600)     # wait 10 minutes

        if (wait is not None and
            submission_count == wait[1] and
            job_count <= jobs_requested):
               print(wait[1] + ' jobs submitted. Pausing for ' +
                     wait[0] + ' minutes.')
               submission_count = 0
               time.sleep(wait[0] * 60)

        if job_count <= jobs_requested:
            submission = ggdc_submit(url, email, blastVariant, qfile, rfile)
            job_count += 1
            submission_count += 1
            if 'Dear User,' in submission:
                print(submission)
                print('Successfully submitted job ' + str(job_count) + '.')
                if wait is not None:
                    print('This is job ' + str(submission_count) +
                          ' of this submission set.')
            else:
                print(submission)
                print('Job ' + str(job_count) +
                      ' failed. Skipping to the next job.')
            time.sleep(2)

        else:
            print(('Error with GGDC job ' + job_count +
                   '. Skipping to the next job.'))

def main(args):

    url = 'http://ggdc.dsmz.de/ggdc.php'

    email = args.email
    blastVariant = args.blastVariant
    queryfile = args.queryfile
    reffile = args.reffile
    samplefile = args.samplefile

    bruteforce = args.bruteforce
    wait = args.wait

    maxrefs = 75

    if args.queryfile and args.reffile is not None:
        pairs_dict = build_pairs_rq(queryfile,reffile)
    elif args.queryfile is not None and args.reffile is None:
        sys.exit(('You must submit a query file and reference file together. '
                  'Exiting.'))
    elif args.samplefile is not None:
        pairs_dict = build_pairs_all(samplefile)
    else:
        sys.exit('Error with data file specification on the command line. Exiting.')

    submissions_dir = os.path.join(get_script_path(), 'submissions')
    if not os.path.exists(submissions_dir): os.makedirs(submissions_dir)
    files_dict = write_submission_files(pairs_dict,
                                        submissions_dir,
                                        maxrefs = 75)

    submit_ggdc_jobs(url, email, blastVariant, files_dict, bruteforce, wait)

if __name__ == "__main__":
    main(args)
