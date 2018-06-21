#!/usr/bin/env python

# An updated and generalized version of scripts originally developed
# by Katya McGough at American Type Culture Collection

# Automatically checks server load on the Genome-to-Genome Distance Calculator
# (GGDC) website: https://ggdc.dsmz.de/ggdc.php

import sys
import os
import mechanize
from bs4 import BeautifulSoup

url = 'http://ggdc.dsmz.de/ggdc.php'

# browse to GGDC website
br = mechanize.Browser()
page = br.open(url)

# check server load
page_html = BeautifulSoup(page, 'html.parser')
status_html = page_html.find('div', {'class': 'progress'})
status = status_html.get_text()
status_message = 'Current GGDC server load:' + status
print(status_message)
