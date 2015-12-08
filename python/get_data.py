#! /usr/bin/env python3
# encoding utf-8

import sys
import functools
import itertools
from collections import defaultdict
import pickle
import argparse
import random

from Bio import Entrez
import socket
from pyphy import *

import atexit
from time import clock, sleep

def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        functools.reduce(lambda ll,b : divmod(ll[0],b) + ll[1:],
            [(t*1000,),1000,60,60])

line = "="*40
def log(s, elapsed=None):
    print(line)
    if elapsed:
        print("Execution time: %s" % elapsed)
    print(line)

def endlog():
    end = clock()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))

def now():
    return secondsToStr(clock())

start = clock()
atexit.register(endlog)

import logging
if "logger" not in globals():
    logger = logging.getLogger('kmerize')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(filename)s - %(message)s',"%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)


alphabet = ['A', 'C', 'G', 'T']
joker = '*'
pattern = '1101'

categories = {0:'Archaea', 1:'Bacteria', 2:'Eukaryota', 3:'Viruses'}
chunks = 50 * len(categories)
read_size=1000

Entrez.email = 'marie.gasparoux@u-bordeaux.fr'

host = 'localhost'
port = 1234
sleep_delay = 5
save_sign = 'save'

pickle_file = "species_gis_by_cat.pickle"

def space_kmer(kmer, pattern):
    return "".join(letter if bit=='1'
                   else joker
                   for letter, bit in zip(kmer, pattern))

def kmerize(read, pattern, read_size):
    d = defaultdict(int)
    lp = len(pattern)
    max_pos = len(read)-lp+1
    for i in range(max_pos):
        d[space_kmer(read[i:i+lp], pattern)] += 1
    coef = len(read)/read_size
    return {k:v*coef for k, v in d.items()}

def print_vw_features(g_label, sk_dict, g_cat=None):
    if g_cat != None:
        feats = "%d '%d |" % (g_cat, g_label)
    else:
        feats = "'%d |" % g_label
    for k, v in sk_dict.items():
        feats += " %s:%d" % (k, v)
    return feats + "\n"

def group(iterable, n):
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args)

def get_category_gis(cat):
    logger.debug("- Category: %s" % cat)
    logger.debug("- Getting all taxids...")
    taxids = getAllSonsByTaxid(getTaxidByName(cat)[0])
    logger.debug("- Converting species taxids to gis...")
    return [gi_list for gi_list in (getGiByTaxid(taxid) for taxid in taxids if getRankByTaxid(taxid) == 'species') if gi_list]

def get_NCBI_sequence(gis):
    ids_list = ",".join([str(gi) for gi in gis if gi])
    handle = Entrez.efetch('nuccore', id=ids_list, rettype='fasta', retmode='xml')
    records = Entrez.parse(handle)
    gi_seq = [(record['TSeq_gi'], record['TSeq_sequence']) for record in records]
    handle.close()
    return gi_seq

def main(argv=None):
    global args
    parser = argparse.ArgumentParser(description='Feed Vowpal Wabbit with contigs from NCBI sequences')
    parser.add_argument('-p', dest='pickle',
                        help='Read GI list from pickle file instead of creating it', action='store_true')
    args = parser.parse_args()

    logger.debug("PREPROCESSING")
    logger.debug("Getting all gis...")
    if not args.pickle:
        gis = [get_category_gis(categories[k]) for k in categories.keys()]
        with open(pickle_file, 'wb') as output:
            pickle.dump(gis, output)
    else:
        with open(pickle_file, 'rb') as input:
            gis = pickle.load(input)

    logger.debug("TRAINING")
    logger.debug("Selecting some gis...")

    random.seed()

    batch = 0
    
    train_gis = [random.sample(gis[k], 500) for k in categories.keys()]
    grped_alt_cat_gis = group(itertools.chain.from_iterable(zip(train_gis[0], train_gis[1], train_gis[2], train_gis[3])), chunks)
    
    for gi_list in grped_alt_cat_gis:
        batch += 1
        logger.debug("Batch %d" % batch)
        logger.debug("Downloading sequences...")
        gi2seq = get_NCBI_sequence(gi_list)
        logger.debug("Kmerizing and sending to vw...")
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((host, port))
        for i in range(len(gi2seq)):
            (gi, seq) = gi2seq[i]
            features = kmerize(seq[:read_size], pattern, read_size)
            vw_feat = print_vw_features(int(gi), features, i%4 + 1)
            s.sendall(vw_feat.encode())
            
        s.sendall(save_sign.encode())
        s.shutdown(socket.SHUT_WR)
        logger.debug("Sleeping (%d s)" % sleep_delay)
        sleep(sleep_delay)

    return 0

if __name__ == '__main__':
    sys.exit(main())
