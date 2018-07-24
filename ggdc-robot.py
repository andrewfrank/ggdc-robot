#!/usr/bin/env python

__version__ = '0.0.4'

# Automatically submit jobs to the Genome-to-Genome Distance Calculator (GGDC)
# website: https://ggdc.dsmz.de/ggdc.php, using GGDC v2.1.

# DEPENDENCIES

import argparse
import os
import sys
import itertools
import numpy
import time
import requests

# COMMAND LINE ARGUMENTS

description = ('Automatically submit jobs to the Genome-to-Genome'
               'Distance Calculator (GGDC) website:'
               'https://ggdc.dsmz.de/ggdc.php. Command line arguments must '
               'include either BOTH --queryfile AND --reffile, OR ONLY'
               '--samplefile.')

parser = argparse.ArgumentParser(description = description)
ifiles = parser.add_mutually_exclusive_group(required = True)
submit = parser.add_mutually_exclusive_group(required = True)

# Required file arguments
ifiles.add_argument('--samplefile','-s',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for all sample '
                    'sequences OR the full path to a fasta file containing '
                    'sample sequences; identity for the entire list is parsed '
                    'from the first entry; NOTICE: THIS OPTION IS '
                    'POTENTIALLY EXTREMELY COMPUTATIONALLY INTENSIVE.'),
                    metavar = 'FILE')
ifiles.add_argument('--queryfile','-q',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a query '
                    'sequence OR the full path to a fasta file containing '
                    'query sequence(s); identity for the entire list is '
                    'parsed from the first entry.'),
                    metavar = 'FILE')
parser.add_argument('--reffile','-r',
                    help = ('the full path to a text file where each new line '
                    'contains EITHER the NCBI accession number for a reference '
                    'sequence OR the full path to a fasta file containing '
                    'reference sequence(s); identity for the entire list is '
                    'parsed from the first entry.'),
                    metavar = 'FILE',
                    required = '--queryfile' in sys.argv)

# Required arguments
parser.add_argument('--email','-e',
                    help = 'the email address where GGDC will send results.',
                    metavar = 'EMAIL ADDRESS',
                    required = True)
submit.add_argument('--slotusage', '-u',
                    help = ('enable slot usage waiting mode; this mode forces '
                    'ggdc-robot to pause for 10 minutes before attempting to '
                    'submit another job when GGDC servers reach this '
                    'specificied capacity. E.g. --slotusage 50 will prompt '
                    'ggdc-robot to wait 10 minutes when GGDC server slot usage '
                    'reaches 50 percent.'),
                    metavar = 'PERCENT')
submit.add_argument('--bruteforce', '-f',
                    help = ('enable brute force mode; this mode forces '
                    'ggdc-robot to submit jobs to the GGDC server even when the'
                    'server load is at 100 percent; ATTENTION: jobs may fail due to GGDC job queue limits.'),
                    action = 'store_true')

# Optional arguments
parser.add_argument('--blastVariant','-b',
                    help = ('the alignment tool used to determine matches '
                    'between query and reference genomes; GGDC recommends '
                    'GBDP2_BLASTPLUS.'),
                    choices = ['GBDP2_BLASTPLUS','GBDP2_BLAT','GBDP2_BLASTZ',
                    'GBDP2_WU-BLAST','GBDP2_MUMMER'],
                    default = 'GBDP2_BLASTPLUS')
parser.add_argument('--timedwait', '-t',
                    help = ('enable timed waiting mode; this mode forces '
                    'ggdc-robot to wait the X minutes every Y jobs '
                    'submitted. E.g. --wait 25 6 will force ggdc-robot to '
                    'wait 25 minutes between every set of 6 jobs submitted. '
                    'Can be combined with --usagewait if desired.'),
                    metavar = ('MINUTES','JOBS'),
                    type = int,
                    nargs = 2)

# Help arguments
# parser.add_argument('--verbose','-v',
#                     help = ('outputs text based checkpoints and interactive '
#                     'submission check.'),
#                     action = 'store_true')
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
    page_request = requests.get(url)

    # check server load
    page_html = page_request.text
    status_loc = page_html.find('aria-valuenow=') + 15
    status = page_html[status_loc]
    return(status)

def ggdc_submit(url, email, blastVariant, queryfile, reffile):

    # begin filling out GGDC form
    form = [('blastVariant',(None,blastVariant)),
            ("targetName",(None,"")),
            ("targetGenome", (None, "")),
            ("refGenbank",(None,"")),
            ("multipleRefGenomes[]",(None, "")),
            ('email',(None,email)),
            ('singlebutton',(None,""))]

    # fill in query form from queryfile
    queryformat = check_submission_format(queryfile)
    with open(queryfile) as infile: qlines = infile.read().splitlines()
    qline = qlines[0]

    if queryformat == 'accession':
        form_value = ("targetName",(None,qline))
        form[1] = form_value
    elif queryformat == 'filepath':
        form_value = ("targetGenome",
                      (qline,open(qline,"rb"),"application/octet-stream"))
        form[2] = form_value
    else: sys.exit('Error submitting query' + queryformat + '. Exiting.')

    # fill in ref form from reffile
    refformat = check_submission_format(reffile)
    with open(reffile) as infile: rlines = infile.read().splitlines()

    if refformat == 'accession':
        control = 'refGenbank'
        sep = '\r\n'
        values = sep.join(rlines)
        form_value = ("refGenbank",(None,values))
        form[3] = form_value
    elif refformat == 'filepath':
        mrg_form_values = []
        for rline in rlines:
            form_value = ("multipleRefGenomes[]",
                          (rline, open(rline,"rb"),"application/octet-stream"))
            mrg_form_values.append(form_value)
        form[4] = mrg_form_values.pop(-1)
        for v in mrg_form_values:
            form.insert(4, mrg_form_values.pop(-1))
    else: sys.exit('Unable to submit reference' + refformat + '. Exiting.')

    # submit GGDC job
    submission = requests.post(url, files = form)

    # get GGDC job response
    submission_html = submission.content.decode()
    if 'job with ID' in submission_html: response = 'Succeeded'
    else: response = 'Failed'
    return(response)

# iteratively submits each qfile-rfile pair to GGDC using ggdc-crawler.py;
# currently pauses for 25 minutes every 6th submission
def ggdc_submission_controller(status_url, submit_url, email, blastVariant,
                               files_dict, bruteforce, wait, slotusage):

    submission_count = 0
    jobs_requested = len(files_dict)
    print('Jobs requested = ' + str(jobs_requested))
    print('---------------------------------------------------')

    for job_count, (qfile, rfile) in enumerate(files_dict.items(), 1):
        job_counter = str(job_count) + '/' + str(jobs_requested)

        status = get_ggdc_status(status_url)
        status_message = 'Current GGDC server load:' + str(status) + '%'
        print(status_message)

        if bruteforce is False:
           if status >= slotusage:
               print( slotusage +
                     (' percent of GGDC server slots are currently used, '
                      'waiting 10 minutes before attempting additional '
                      'submissions.'))
               print('Waiting to submit job ' + job_counter + '.')
               if wait is not None:
                   print('This is job ' + str(submission_count) +
                         ' of this submission set.')
               print('---------------------------------------------------')
               time.sleep(600)     # wait 10 minutes

        if (wait is not None and
            submission_count == wait[1] and
            job_count <= jobs_requested):
               print(wait[1] + ' jobs submitted. Pausing for ' +
                     wait[0] + ' minutes.')
               submission_count = 0
               time.sleep(wait[0] * 60)

        if job_count <= jobs_requested:
            submission = ggdc_submit(submit_url, email, blastVariant,
                                     qfile, rfile)
            job_count += 1
            submission_count += 1
            print('GGDC server submission attempt:')
            if 'Succeeded' in submission:
                print(submission)
                print('Successfully submitted job ' + job_counter + '.')
                if wait is not None:
                    print('This is job ' + str(submission_count) +
                          ' of this submission set.')
                print('---------------------------------------------------')

            else:
                print(submission)
                print('Job ' + job_counter +
                      ' failed. Skipping to the next job.')
                print('---------------------------------------------------')
            time.sleep(2)

        else:
            print(('Error with GGDC job ' + job_counter +
                   '. Skipping to the next job.'))
            print('---------------------------------------------------')

def main(args):

    status_url = 'http://ggdc.dsmz.de/ggdc.php'
    submit_url = 'http://ggdc.dsmz.de/submit_ggdc_job.php'

    email = args.email
    blastVariant = args.blastVariant
    queryfile = args.queryfile
    reffile = args.reffile
    samplefile = args.samplefile

    slotusage = args.slotusage
    bruteforce = args.bruteforce
    wait = args.timedwait

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

    ggdc_submission_controller(status_url, submit_url, email, blastVariant,
                               files_dict, bruteforce, wait, slotusage)
    sys.exit('Script complete. Exiting.')

if __name__ == "__main__":
    main(args)
