# encoding utf-8

import sys
import pickle
from pyphy import *

categories = {0:'Archaea', 1:'Bacteria', 2:'Eukaryota', 3:'Viruses'}

def main(argv=None):
    print("PREPROCESS")
    print("Getting taxids...")

    taxids = [[taxid for taxid
               in getAllSonsByTaxid(getTaxidByName(cat)[0])
               if getRankByTaxid(taxid) == 'species' ]
              for cat in categories.values()]

    print([len(l) for l in taxids])

    with open("taxids.pickle", 'wb') as output:
        pickle.dump(taxids, output)
    return 0

if __name__ == '__main__':
    sys.exit(main())

