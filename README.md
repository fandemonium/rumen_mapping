# shallow sequenced metagenome analysis
## the project has a emphasis on rumen but could be applied to various environment. Most importantly, finding the most relevant genome data base. 
## If all fails, there is always the RefSeq and blast ...

### reference data set: 
+ from paper: https://www.nature.com/articles/s41467-018-03317-6#Sec17
 
  `Assembly of 913 microbial genomes from metagenomic sequencing of the cow rumen`

+ assembled genome and proteome files can be downloaded via:
  ```
  cd /PATH/TO/WHERE/EVERYTHING/IS
  mkdir rumen && cd rumen
  wget http://datashare.is.ed.ac.uk/download/DS_10283_2772.zip
  ```
  + NOTE: need to do multiple unzip and untar.

+ the datashare site above only contains wgs contigs, and protein sequences. The identified and binned rumen unculatured genomes (RUG) are in ENA (as well as the raw sequences etc...). 
  ```
  # ENA study accession number: PRJEB21624
  # use ENA protal api generater: will print a table of what the study contains
  curl -X GET --header 'Accept: text/plain' 'https://www.ebi.ac.uk/ena/portal/api/links/study?accession=PRJEB21624' 
  ####result_id	description	entry_cnt
  ####assembly	Assembly	913
  ####read_experiment	Experiment	45
  ####read_run	Read	45
  ####wgs_set	Genome assembly contig set	913
  # then set the results for RUG IDs (may need to surpress cert error `-k`): wgs_acc RUG id
  curl -k -X GET --header 'Accept: text/plain' 'https://www.ebi.ac.uk/ena/portal/api/links/study?accession=PRJEB21624&result=wgs_set' > PRJEB21624_wgs_set.txt
  # the first column is the RUG wgs id (almost... needs to get rid of the trailing 0's)
  cd rumen
  mkdir RUG_genomes & cd RUG_genomes
  while read line; do curl -O http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/om/$line.fasta.gz; done < <(tail -n +2 ../PRJEB21624_wgs_set.txt | cut -f 1 | sed 's/000000$//g')
  # should 913 of them. and the headers contain RMG contig/protein ids. 
  ```
  
  
### analysis prcedures:

+ use bwa.
  + create idex for the reference genome (need to concat all and checked. no duplicated headers):
    ```
    cd rumen/genomes
    cat *.fa > all_rmgs.fa
    bwa index -p ../bwa_index/rmg_genomes all_rmgs.fa
    ```
    
    NOTE: don't use pipe to combine cat and bwa... it runs into issues.
  
  + location of the shallow sequenced metagenomes
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd 181214_fastqs
    ```
    
  + quality filter fastqs: (maxee 0.5)
    ```
    mkdir filtered_fq && cd fastq  
    for i in *.gz; do vsearch --threads 4 --fastq_filter $i --fastq_maxee 0.5 --fastq_maxns 0 --fastqout ../filtered_fq/${i//.fastq.gz/.maxee0.5.fq}; done
    ```
  
  + mapping
    ```
    mkdir bwa_bams
    cd filtered_fq
    for i in *.fq; do bwa mem -t 4 /mnt/scratch/yangfan1/rumen/bwa_index/rmg_genomes $i | samtools sort -@4 -o ../bwa_bams/${i//.fq/.sorted.bam} -; done
    ```

  + get mapped reads only:
    ```
    mkdir ../bwa_bam
    for i in *.bam; do samtools view -b -F 4 -q 10 $i > ../mapped_bam/${i//sorted/mapped}; done
    ```
    
  + OPTIONAL: mpileup (reference fa):
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd 181214_fastqs
    cd mapped_bam
    mkdir ../ref_mpileup
    for i in *.bam; do samtools mpileup -f ../../rumen/genomes/all_rmgs.fa $i > ../noref_mpileup/${i//.bam/.ref.txt}; done
    ```
  
  + get reference header:
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd rumen/genomes
    grep ">" RMG_*.fa > ../RMG_genome_to_contigs.txt
    cd ../proteomes
    grep ">" RMG_*.faa > ../RMG_genome_to_proteome.txt
    ```

  + create bed file from the proteome txt file:
    ```
    paste -d "\t" <(cut -d ">" -f 2 RMG_genome_to_proteome.txt | cut -d " " -f 1 | rev | cut -d "_" -f 2- | rev) <(cut -d " " -f 3 RMG_genome_to_proteome.txt) <(cut -d " " -f 5 RMG_genome_to_proteome.txt) <(cut -d " " -f 1 RMG_genome_to_proteome.txt) <(cut -d " " -f 7 RMG_genome_to_proteome.txt) > RMG_proteome.bed
    sed -i 's/.faa//g' RMG_proteome.bed 
    sed -i 's/>//g' RMG_proteome.bed
    rev RMG_proteome.bed | cut -f 1 | rev | sed 's/\-1/\-/g' | sed 's/1/\+/g' > temp.txt
    paste -d "\t" <(cut -f 1-4 RMG_proteome.bed) <(rev RMG_genome_to_proteome.txt | cut -d "=" -f 1 | rev) <(cut -f 5 RMG_proteome.bed) > temp.bed
    mv rumen/temp.bed rumen/RMG_proteome.bed
    ```
    
      bed file in a format like this (need to be 6 columns): id start end name some_value(eg. gc_content) strand
      ```
      k87_58769312    514     1122    RMG_1025:k87_58769312_1 0.268   +
      k87_58769312    1307    2620    RMG_1025:k87_58769312_2 0.215   +
      k87_58269671    207     698     RMG_1025:k87_58269671_1 0.250   -
      ```
  
  + bedtools to find intersect and coverage:
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd 181214_fastqs
    mkdir intersects
    cd mapped_bam
    # find intersect:
    for i in *.bam; do bedtools intersect -a ../../rumen/RMG_proteome.bed -b $i -bed -wa -wb -s > ../intersects/${i//mapped.bam/intersect.bed}; done
    # find coverage
    mkdir ../coverages
    for i in *.bam; do bedtools coverage -a ../../rumen/RMG_proteome.bed -b $i -s | grep -vw "0.0000000" > ../coverages/${i//mapped.bam/coverage.txt}; done
    # coverage hist may be useful too.
    mkdir ../coverages_hist
    for i in *.bam; do bedtools coverage -a ../../rumen/RMG_proteome.bed -b $i -s -hist | grep -vw "1.0000000" > ../coverages_hist/${i//mapped.bam/coverage_hist.txt}; done
    ```
    
    
+ use blastx (diamond)
 
  + get nr:
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    mkdir nr_20190414 && cd nr_20190414/
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.gz
    for i in *.gz; do tar -zxvf $i; done
    ```
