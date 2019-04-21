import sys
import gzip
import re
from Bio import SeqIO

with gzip.open(sys.argv[1]) as f:
	for record in SeqIO.parse(f, "fasta"):
		bwa_id = record.id
		rug_id = "RUG" + re.split("RUG", record.description)[1].split(" ")[0]
		contig = record.description.strip().split(" ")[-1]
		print bwa_id + "\t" + "1" + "\t" + str(len(record.seq)) + "\t" + rug_id + ":" + contig
