import sys
import re

for lines in open(sys.argv[1], 'rU'):
	line = lines.strip()
	bwa_name = re.split(">| ", line)[1]
	rug_id = "RUG" + re.split("RUG", line)[1].split(" ")[0]
	contig = line.split(" ")[-1]
	print bwa_name + "\t" + rug_id + "\t" + contig
