# shallow sequenced metagenome analysis
## the project has a emphasis on rumen but could be applied to various environment. Most importantly, finding the most relevant genome data base. 
## If all fails, there is always the RefSeq...

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

### analysis prcedures:

+ use bwa.
  + create idex for the reference genome:
    ```
    bwa index -p ../bwa_index/rmg_genomes *.fa
    ```
  
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
    for i in *.fq; do echo "bwa mem -t 4 /mnt/scratch/yangfan1/rumen/bwa_index/rmg_genomes $i | samtools sort -@4 -o ../bwa_bams/${i//.faq/.sorted.bam} -"; done > ../bwa_mapping.sh
    bash ../bwa_mapping.sh
    ```

  + get mapped reads only:
    ```
    mkdir ../mapped_bam
    for i in *.bam; do samtools view -b -F 4 $i > ../mapped_bam/${i//sorted/mapped}.bam; done
    ```
    
  + mpileup (without reference fa, since SNP is not important here):
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd mapped_bam
    mkdir ../noref_mpileup
    for i in *.bam; do samtools mpileup $i > ../noref_mpileup/${i//.bam.bam/.min3.txt}; done
    ```
    
+ use blastx (diamond)
 
  + get nr:
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    mkdir nr_20190414 && cd nr_20190414/
    wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.gz
    for i in *.gz; do tar -zxvf $i; done
    ```
