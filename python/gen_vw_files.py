from Bio import SeqIO
from Bio.Alphabet import generic_dna
from get_data import *

categories = {0:"archea", 1:"bact", 2:"euk", 3:"viruses"}
output = open("15000contigs_gi_pos_cat", 'w')
pattern = '11011'

for i in categories.keys():
    with open(categories[i] + "_" + pattern, 'w') as output:
        with open(categories[i] + "_15000contigs.bfast.fastq", 'r') as file:
            records = SeqIO.parse(file, "fastq")
            id = 0
            second_read = False
            for r in records:
                if not second_read:
                    seq = str(r.seq)
                    vect = kmerize(seq, pattern, 1000)
                    feat = print_vw_features(id, vect, i)
                    output.write(feat)
                    id += 1
                second_read = not second_read
