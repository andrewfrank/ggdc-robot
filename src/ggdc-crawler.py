#!/usr/bin/env python

# An updated and generalized version of scripts originally developed
# by Katya McGough at American Type Culture Collection

# Automatically submit jobs to the Genome-to-Genome Distance Calculator (GGDC)
# website: https://ggdc.dsmz.de/ggdc.php, using GGDC v2.1

# todo
# - handle ggdc server errors properly
# - implement submission of multiple reference genome fasta files

# notes
# - currently using this search criteria on genbank to find genome accessions:
# "Escherichia"[Organism] AND (complete[Properties] or "wgs master"[Properties]) 

# DEPENDENCIES

import sys
import os
import mechanize
from bs4 import BeautifulSoup

# FUNCTIONS

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

# SCRIPT

url = 'http://ggdc.dsmz.de/ggdc.php'
email = sys.argv[1]
blastVariant = [sys.argv[2]]
queryfile = sys.argv[3]
reffile = sys.argv[4]

br = mechanize.Browser()
br.open(url)
br.select_form('Form')                          # GGDC form name is just "Form"

br.form.set_value(email, 'email')               # fill in email accession form
br.form.set_value(blastVariant, 'blastVariant') # fill in GGDC BLAST method form

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
# change this to handle error messages from ggdc as well
message_html = submission_html.find('div', {'class': 'alert alert-success'})
message = message_html.get_text()
print(message)
