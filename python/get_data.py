# encoding utf-8

import logging, timing, time
logger = logging.getLogger('metaviro_logger')

import sys, os
import itertools
from collections import defaultdict
import argparse
import random
import pickle

import ncbi_data, vw_daemon
from pyphy import *

categories = {0:'Archaea', 1:'Bacteria', 2:'Eukaryota', 3:'Viruses'}
ref_read_size=1000

sleep_delay = 1


def space_kmer(kmer, pattern):
    return "".join(letter if bit=='1'
                   else '*'
                   for letter, bit in zip(kmer, pattern))


def kmerize(read, pattern, ref_read_size):
    d = defaultdict(float)
    lp = len(pattern)
    max_pos = len(read)-lp+1
    for i in range(max_pos):
        d[space_kmer(read[i:i+lp], pattern)] += 1
    coef = len(read)/ref_read_size
    return {k:v*coef for k, v in d.items()}


def random_sampling(taxids):
    train_taxids = random.sample(taxids, 200)
    train_gis = [random.choice(gi_list) for gi_list
                    in [getGiByTaxid(taxid) for taxid in train_taxids]
                    if gi_list != []]
    grp_size = 20
    final_train_gis = []
    grped_train_gis = [train_gis[i:i+grp_size] for i in range(0,len(train_gis),grp_size)]
    for gi_list in grped_train_gis:
        try:
            ids_list = ",".join([str(gi) for gi in gi_list])
            records = ncbi_data.get(ncbi_data.get_gis_pos, ids_list)
            i = 0
            for record in records:
                rlen = record['Length']
                if rlen <= ref_read_size:
                    random_pos = 0
                    read_size = rlen
                else:
                    random_pos = random.randint(0, rlen - ref_read_size)
                    read_size = ref_read_size
                final_train_gis.append([gi_list[i], random_pos, read_size])
                i += 1
        except RuntimeError as e:
            logger.error("GI not found in database: %s", str(e))
    return final_train_gis


def main(argv=None):
    parser = argparse.ArgumentParser(description='Download reference sequences from NCBI datasets and feed the vw classifier(s)')
    parser.add_argument('--pattern_port', dest='pattern_port_file', type=str, required=True, 
                        help='Tabular file listing for each daemon: initial_model, pattern, port and destination repertory. See config/model_pattern_port_dest.txt for an example')
    parser.add_argument('--debug', dest='debug', action='store_true', 
                        help='For detailed output')
    args = parser.parse_args()

    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logger.setLevel(log_level)
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    formatter = logging.Formatter("%(asctime)s (%(levelname)s) - %(message)s",
                                  "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    logger.info("PREPROCESS")

    with open(args.pattern_port_file, 'r') as input:
        pattern_port_dest_list = [[info for info in line.split("\n")[0].split("\t")[1:] if info]
                                  for line in input]
    logger.debug("Parsing pattern_port file for [pattern, port, dest]")
    logger.debug(pattern_port_dest_list)
    pattern_port_dict = defaultdict(list)
    port_dest_dict = {}
    for pattern_port in pattern_port_dest_list:
        pattern_port_dict[pattern_port[0]].append(pattern_port[1])
        port_dest_dict[pattern_port[1]] = pattern_port[2]

    with open("pattern_port_dict.pickle", 'wb') as output:
        pickle.dump(pattern_port_dict, output)
    with open("port_dest_dict.pickle", 'wb') as output:
        pickle.dump(port_dest_dict, output)

    random.seed()
    train_gis = [[]]*4
    cpt = 0

    logger.info("Getting taxids...")
    with open("taxids.pickle", 'rb') as input:
        taxids = pickle.load(input)

    logger.info("TRAINING")
    while True:
        for k in categories.keys():
            train_gis[k] = random_sampling(taxids[k])
                
        min_size = min([len(l) for l in train_gis])
        cpt += min_size
        logger.debug([len(l) for l in train_gis])
        train_gis = [l[:min_size] for l in train_gis]
        logger.debug(train_gis)

        for i in range(min_size):
            for cat in categories.keys():
                (gi, pos, read_size) = train_gis[cat][i]
                logger.debug((cat, gi, pos, read_size))
                seq = ncbi_data.get(ncbi_data.get_gis_slice_seq, gi, pos, pos+read_size-1)
                logger.debug(seq)
                for pattern in pattern_port_dict.keys():
                    logger.debug(pattern)
                    features = kmerize(seq, pattern, ref_read_size)
                    vw_feat = vw_daemon.print_vw_features(int(gi), features, cat + 1)
                    logger.debug(pattern_port_dict[pattern])
                    logger.debug(vw_feat)
                    for port in pattern_port_dict[pattern]:
                        vw_daemon.send_to_daemon(port, vw_feat)
            time.sleep(sleep_delay)
        logger.info(cpt)

    return 0

if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        try:
            print()
            logger.info("Interrupted by user")
            sys.exit(0)
        except SystemExit:
            os._exit(0)
