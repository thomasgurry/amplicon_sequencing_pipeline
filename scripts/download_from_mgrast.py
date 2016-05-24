#!/usr/bin/env python
# This script calls the MG-RAST API and attempts to download the raw data file for a
# metagenome whose MG-RAST ID is specified.
import sys
import os
import re
from optparse import OptionParser
import urllib2
# MG-RAST API url
API_URL = "http://api.metagenomics.anl.gov/1"
def retrieveMGRbyaccession(accessionno, key):
    '''Retrieve raw data from MG-RAST API using curl and dump result into file named <accession>.gz'''
    try:
        a = re.search(r"(\d\d\d\d\d\d\d\.\d)$", accessionno).group(1)
    except IndexError:
        sys.exit("Don't recognize accession number format %s" % accessionno)
    if key == " ":
        s1 = "%s/download/mgm%s?file=050.1" % (API_URL, a)
    else:
        s1 = "%s/download/mgm%s?file=050.1&auth=%s" % (API_URL,a,key)
    sys.stderr.write("Retrieving %s\n" % s1)
    try:
        opener = urllib2.urlopen(s1)
        open("%s.gz" % a, "w").write(opener.read())
    except urllib2.HTTPError, e:
        print "Error with HTTP request: %d %s\n%s" % (e.code, e.reason, e.read())
        sys.exit(255)


if __name__ == '__main__':
    usage = '''usage: download_metagenome_sequences.py <accession number>
example: download_metagenome_sequences.py MGR4440613.3'''
    parser = OptionParser(usage)
    (opts, args) = parser.parse_args()
# Assign the value of key from the OS environment
    try:
        key = os.environ["MGRKEY"]
        #key='ccjdmXCmxePvPHHFauF7yhDuY'
#        key='eFr8ksKuehV6eeFfCdV6PC7LA'
    except KeyError:
        key = ""
    if len(key) == 25:
        sys.stderr.write("Using MGR webkey %s\n" % key)
    else:
        sys.stderr.write("Warning: MGR webkey not defined\n")
# test for correct number of arguments
    try:
        accession = args[0]
    except IndexError:
        parser.error("accession is a required parameter\n")

# call retrieveMGRbyaccession
    retrieveMGRbyaccession(accession, key)
