#!/usr/bin/env python3

from collections import defaultdict
import argparse, gzip, sys

def kmerize(read, pattern, ref_read_size):
	
	d = defaultdict(float)
	lp = len(pattern)
	max_pos = len(read)-lp+1
	for i in range(max_pos):
		d[space_kmer(read[i:i+lp], pattern)] += 1
	coef = len(read)/ref_read_size
	return {k:v*coef for k, v in d.items()}


def print_vw_features(g_label, sk_dict, g_cat=None):
	
	if g_cat != None:
		feats = "%d '%d |" % (g_cat, g_label)
	else:
		feats = "'%d |" % g_label
	for k, v in sorted(sk_dict.items()):
		feats += " %s:%.2f" % (k, v)
	return feats + "\n"

def parse_function():
	print("parsing arguments")
	parser=argparse.ArgumentParser(description="Calculates k-mer frequencies from a given fasta file(s).")
	parser.add_argument('-f', '--file', dest="infile", type =str, help="fasta file for the sequences (left empty to read from stdin)")
	parser.add_argument("-k", "--k-mer", dest="kmer", type=int, help="the pattern of the given k-mer like 111")
	parser.add_argument("-c", "--category", dest="category", type=int, help="The category of the given sequences")
	parser.add_argument("-l", "--label", dest="label", type=str, help="Human readable label information for the category")
	parser.add_argument("-o", "--out-file",dest="outfile",type=str,help="name of the resulting file")
	args=parser.parse_args()

	print("parsed")

	if len(sys.argv)==1:
		print("no args")
		# display help message when no args are passed.
		parser.print_help()
	sys.exit(1)

	return(args)           


def main():

	print("main")
	args=parse_function()
	infile=args.infile
	kmer=str(args.kmer)
	category=args.category
	label=args.label
	outfile=args.outfile
	
	print("trying to open")
	print(infile)
	infile="set.fasta"
	print(infile)
	print("amk")


	try:
		extension=infile[-2:]
		if extension=="gz":
			print(infile)
			seq_file=gzip.open(infile,"r")
		else:
			seq_file=open(infile, 'r')
			print(seq_file)
	except:
		seq_file=open("/dev/stdin","r")
		print(seq_file)

	line=seq_file.readline()
	print(line)
	if line[0]!=">":
		print("Not a fasta file")
		exit()
	output=open(outfile,'w')
	print("file opened")

	while line!='':

		if line[0]==">":
			seq=''
			id=line.split()[0][1:]
			line=seq_file.readline().rstrip("\n")
		if line[0]!=">":
			seq=seq+line
			line=seq_file.readline().rstrip("\n")

		if line=="":
			features=kmerize(seq,kmer,150)
			print(print_vw_features(category,features,label), file=output)
			break
		if line[0]==">":
			features=kmerize(seq,kmer,150)
			print(print_vw_features(category,features,label), file=output)

	output.close()	

if __name__ == '__main__':                  #the main function
	print("yayyy!")
	main()
