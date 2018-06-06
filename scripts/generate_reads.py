#!/usr/bin/env python


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
import re, sys, csv, StringIO, random, decimal, argparse

parser = argparse.ArgumentParser(description='Generate uniform-length single or paired-end metagenomic reads.')
parser.add_argument('-r', metavar='<reference_fasta>', dest="ref", help="Multi-FASTA file containing genomic sequences from which reads will be sampled.")
parser.add_argument('-a', metavar='<abundance_file>', dest="abund", help="Tab-delimited abundance file with an abundance value for each corre- sponding genome sequence in <reference fasta>")
parser.add_argument('-o', metavar='<output_file>', dest="output", help="Name for output file containing simulated uniform-length reads")
parser.add_argument('-t', metavar='<total_reads>', type=int, dest="total", help="The total number of reads to sample from all genomes")
parser.add_argument('-l', metavar='<longest_read>', type=int, dest="length", help="The length, in bp, of the longest possible read to simulate")
parser.add_argument('-i', metavar='<insert_mean_length>', type=int, dest="insert", default="0", help="Average length of insert for paired-end reads.")
parser.add_argument('-s', metavar='<insert_stddev>', type=int, dest="stddev", default="0", help="Standard deviation of insert length for paired-end reads" )
parser.add_argument('-d', '--direction', action='store_true', dest="direction", help="Use this switch to generate reads in both forward and reverse orientations" )
args = parser.parse_args()


# print help if no arguments passed
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit(0)

#Reference metagenome database file (FASTA)
f1 = open(args.ref);

#abundance file (tab-delimited .txt)
f2 = open(args.abund);

total_reads = args.total

max_read_length = args.length

insert_avg = args.insert
insert_stddev = args.stddev

if(insert_avg):
	f4 = open(args.output + '.1.fasta', 'w')
	f5 = open(args.output + '.2.fasta', 'w')
else:
	f4 = open(args.output, 'w')

species=[]
diversity=[]
lengths=[]
freqs=[]

div_file = csv.reader(f2, delimiter='\t')
for row in div_file:
	species.append(row[0][1:])
	diversity.append(decimal.Decimal(row[1]))

uniqueReadCounter = 0

for i in SeqIO.parse(f1, 'fasta') :
	i = i.upper()
	i.seq= Seq(re.sub('[YRWSKMDVHBX]', 'N', str(i.seq)), generic_dna)

	genome_num=0
	while(not(species[genome_num] in i.description)) :
		genome_num+=1
	if(species[genome_num] in i.description) :
		coverage=max(1, int((decimal.Decimal(diversity[genome_num])*total_reads)))

		limit=len(i.seq)
		for j in range(0, coverage) :
			uniqueReadCounter += 1
			rand = random.random()
			rand_length = 0
			numLen = len(lengths)-1
			# making read headers unique
			i.description = i.id + "_read_" + str(uniqueReadCounter)

			if( (insert_avg != 0) & (insert_stddev != 0)):
				cur_insert = int(random.gauss(insert_avg, insert_stddev))
				if(limit > (max_read_length * 2 + cur_insert)):
					start1 = random.randint(0, limit-(2*max_read_length + cur_insert))
					end1 = start1 + max_read_length
					start2 = end1 + cur_insert
					end2 = start2 + max_read_length
				# for contigs upto read length are used and anything smaller are ignored
				elif(limit >= max_read_length):
					start1 = 0
					end1 = max_read_length
					start2 = limit - max_read_length
					end2 = limit
				else:
					continue

				read1 = i.seq[start1:end1]
				read2 = i.seq[end2:start2:-1].reverse_complement()
				# reads from any orienation then randomly make reads from reverse strand
				if(args.direction and random.random() > 0.5):
					read1 = read1.reverse_complement()
					read2 = read2.reverse_complement()
				
				# write reads to files
				f4.write(">%s\n" % i.description)
				f4.write("%s\n" % read1)
				f5.write(">%s\n" % i.description)
				f5.write("%s\n" % read2)

			else:
				if(limit >= max_read_length) :
					start=random.randint(0, limit-max_read_length)
					end=start+max_read_length
				# if contig length less than read length ignore
				else:
					continue
				read = i.seq[start:end]
				#reverse orientation
				if(args.direction and random.random() > 0.5):
					read = read.reverse_complement()

				# write read to file
				f4.write(">%s\n" % i.description)
				f4.write("%s\n" % read)


	if (genome_num >= len(species) ) :
		break;

f1.close()
f2.close()
f4.close()

if(insert_avg):
	f5.close()
