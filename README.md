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

+ use bwa.
  + create idex for the reference genome:
    ```
    bwa index -p ../bwa_index/rmg_genomes *.fa
    ```
  
  + location of the shallow sequenced metagenomes
    ```
    cd /PATH/TO/WHERE/EVERYTHING/IS
    cd 181214_fastqs
    mkdir bwa_bams
    cd fastq
    for i in *.gz; do echo "bwa mem -t 4 /mnt/scratch/yangfan1/rumen/bwa_index/rmg_genomes $i | samtools sort -@4 -o ../bwa_bams/${i//.fastq.gz/.sorted.bam} -"; done > ../bwa_mapping.sh
    bash ../bwa_mapping.sh
    ```


