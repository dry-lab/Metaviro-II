# encoding utf-8
#############################################################################################
#
# Script name : check_NCBI.py
# -----------
# Dev environment : - Ubuntu 14.04 x64
# ---------------   - Python 3.4.3
#
#############################################################################################
#
# Warning : Possible instability due to NCBI servers.
# 
# Use : python3 check_NCBI.py
#
#############################################################################################
#
# History : 2015/11/20 				V0.0 Emmanuel Bouilhol
#
#############################################################################################
#
# Dependencies : - See called packages
#
#############################################################################################
#
# Improvements TBD : 
#
##############################################################################################

## Dependancies
import sys
import functools
import itertools
from collections import defaultdict, Counter
import pickle
import argparse
from Bio import Entrez
import socket
from pyphy import *
import atexit
from time import clock, sleep
import sqlite3
from contextlib import closing
from subprocess import call


## Variables
# SQLITE3 file
db = "../taxdump/ncbi.db"
unknown = -1
no_rank = "no rank"
#Number of days to search
searchLaps = 3
#inititiate tab of GIs
giTab = []
#Hash of taxons
taxHash = {}
#Hash of GIs
giHash = {}
#Entrez parameter
Entrez.email = 'marie.gasparoux@u-bordeaux.fr'

##------------------------------------------------------------------Get total count of new entries
handle = Entrez.esearch(db="nuccore", term="", reldate=searchLaps, datetype="pdat", rettype="count")
record = Entrez.read(handle)
handle.close()
print("New records : ", record["Count"])
count = int(record["Count"])

##------------------------------------------------------------------Get all GIs
i=0
while i <= count:
	handle = Entrez.esearch(db="nuccore", term="", id="gi", reldate=searchLaps, datetype="pdat", retstart=i, retmax=(i+100000))
	record = Entrez.read(handle)
	handle.close()
	for j in range(0, len(record["IdList"])):
		giTab.append(record["IdList"][j])
	i += 100000

##------------------------------------------------------------------Remove existing GIs
curedGiTab = []
for gi in giTab:
	res = getTaxidByGi(gi)
	if res == -1:
		curedGiTab.append(gi)
print("Cured records : ",len(curedGiTab))

##------------------------------------------------------------------Send query to elink by groups of 500 id
print("Retrieving Data from NCBI...")
i=500
while i < len(curedGiTab):
	tempTab = curedGiTab[i-500:i]
	taxHandle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=tempTab)
	taxRecord = Entrez.read(taxHandle)
	taxHandle.close()
	for record in taxRecord:
		if not record["LinkSetDb"]:
			print("NO link to taxonomy")
		else:
			giHash[record['IdList'][0]] = record["LinkSetDb"][0]["Link"][0]["Id"]
			parentTaxId = getParentByTaxid(record["LinkSetDb"][0]["Link"][0]["Id"])
			if parentTaxId == -1:
				taxHash[record["LinkSetDb"][0]["Link"][0]["Id"]] =""
	print (i, "/", len(curedGiTab), end="\r")
	i += 500
print("Number of new GIs : ",len(giHash))
print("Number of new Taxons : ", len(taxHash))

##------------------------------------------------------------------Insert parent taxon in tree table
print("Update taxonomy...")
for taxid in taxHash:
	parentTaxHandle = Entrez.efetch(db="taxonomy", id=taxid, rettype="xml")
	parentTaxRecord = Entrez.read(parentTaxHandle)
	parentTaxHandle.close()
	with sqlite3.connect(db) as conn:
		with closing(conn.cursor()) as cursor:
			cursor.execute("INSERT INTO tree VALUES(\'"+taxid+"\', \'"+parentTaxRecord[0]['ScientificName']+"\', \'"+parentTaxRecord[0]['ParentTaxId']+"\', \'"+parentTaxRecord[0]['Rank']+"\')")

##------------------------------------------------------------------Insert all Gis in gi_taxid table
print("Update GI list...")
for gi, taxid in giHash.items():
	with sqlite3.connect(db) as conn:
		with closing(conn.cursor()) as cursor:
			cursor.execute("INSERT INTO gi_taxid VALUES (\'"+gi+"\', \'"+taxid+"\')")
print("Process complete !")
