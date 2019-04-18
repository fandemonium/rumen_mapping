import sys

with open(sys.argv[1]) as f:
	next(f)
	for lines in f:
		line = lines.strip().split("\t")
		ena_id = line[0].replace("000000", "")
		sub_dir = ena_id[:2].lower()
		print "curl -O http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/" + sub_dir + "/" + ena_id +".fasta.gz"
