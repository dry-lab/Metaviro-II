from Bio import SeqIO
from Bio.Alphabet import generic_dna

categories = {0:"archea", 1:"bact", 2:"euk", 3:"viruses"}

with open("15000contigs_gi_pos_cat", 'w') as output:
    for i in categories.keys():
        with open(categories[i] + "_15000contigs.bfast.fastq", 'r') as file:
            records = SeqIO.parse(file, "fastq")
            second_read = False
            for r in records:
                if not second_read:
                    infos = r.id.split("|")
                    gi = infos[1]
                    read_infos = infos[-1].split("_")
                reverse = int(read_infos[3+second_read])
                if not reverse:
                    pos = read_infos[1+second_read]
                    output.write("%s\t%s\t%d\n" % (gi, pos, i))
                second_read = not second_read
