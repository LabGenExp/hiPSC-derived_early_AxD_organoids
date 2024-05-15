setwd("organoids_experiment1")

#################CREATION WORKING FOLDER###############################################
cat(paste0("mkdir qa;\n",
           "mkdir qa/raw;\n",
           "mkdir qa/trimmomatic;\n",
           "mkdir qa/sortmerna;\n",
           "mkdir qa/fastq_screen;\n",
           "mkdir qa/umitools;\n",
           "mkdir sortmerna;\n",
           "mkdir trimmomatic;\n",
           "mkdir log;\n",
           "mkdir log/dedup;\n",
           "mkdir count;\n",
           "mkdir STAR;\n",
           "mkdir umitools;\n",
           "mkdir count_umi;\n"))

samples<-read.table("samples.csv",header=T,sep=",")
##################QUALITY CONTROL BEFORE ANALYSIS######################################

for(i in 1:nrow(samples)) {
  cat(paste0("zcat ",
             as.vector(samples$LibraryName)[i],"_L001_R1_001.fastq.gz ",
             as.vector(samples$LibraryName)[i],"_L002_R1_001.fastq.gz > ",
             as.vector(samples$SampleName)[i],"_1.fastq;\n",
             "pigz -p 12 ",as.vector(samples$SampleName)[i],"_1.fastq;\n"))
}



cat(paste0("fastqc -t 12 -o qa/raw raw/*"))


################DEFINITION OF VARIABLES################################################

#rRNA databses
ref=paste0("--ref ",
           "/mnt/d/rRNA_database/mtDNA_homo_sapiens.fasta,/mnt/d/rRNA_database/index/mtDNA_homo_sapiens",
           ":/mnt/d/rRNA_database/rfam-5.8s-database-id98.fasta,/mnt/d/rRNA_database/index/rfam-5.8s-database-id98",
           ":/mnt/d/rRNA_database/rfam-5s-database-id98.fasta,/mnt/d/rRNA_database/index/rfam-5s-database-id98",
           ":/mnt/d/rRNA_database/silva-arc-16s-id95.fasta,/mnt/d/rRNA_database/index/silva-arc-16s-id95",
           ":/mnt/d/rRNA_database/silva-arc-23s-id98.fasta,/mnt/d/rRNA_database/index/silva-arc-23s-id98",
           ":/mnt/d/rRNA_database/silva-bac-16s-id90.fasta,/mnt/d/rRNA_database/index/silva-bac-16s-id90",
           ":/mnt/d/rRNA_database/silva-bac-23s-id98.fasta,/mnt/d/rRNA_database/index/silva-bac-23s-id98",
           ":/mnt/d/rRNA_database/silva-euk-18s-id95.fasta,/mnt/d/rRNA_database/index/silva-euk-18s-id95",
           ":/mnt/d/rRNA_database/silva-euk-28s-id98.fasta,/mnt/d/rRNA_database/index/silva-euk-28s-id98")
#HEADCROP
HEADCROP=c(12)
#ADAPTER SEQUENCE
ILLUMINACLIP=c("/mnt/d/Adapters/Lexogen_quantseq.fa:2:30:10")
#GenomeDir
genomeDir=c("/mnt/f/Genomes/Homo_sapiens_GRCh38/")
#GTF file and gene name identificator
GTF=c("/mnt/f/Genomes/Homo_sapiens_GRCh38/Homo_sapiens.GRCh38.87.gtf")
ID=c("gene_id")
#NumberOfThreads
threads=c(10)

#############ANALYSIS SAMPLE AFTER SAMPLE##############################################
#for(i in 71:nrow(samples)) {
for(i in 65:96) {
#for(i in 2:2) {
  #extract umi
  cat(paste0("umi_tools extract",
             " --stdin=raw/",as.vector(samples$SampleName)[i],"_1.fastq.gz",
             " --bc-pattern=\"^(?P<umi_1>.{6})(?P<discard_1>TATA){e<=1}\"",
             " --extract-method=regex",
             " --log=log/",as.vector(samples$SampleName)[i],"_umi-tools.log",
             " --stdout umitools/umi_",as.vector(samples$SampleName)[i],"_1.fastq.gz;\n"))
  
  
  #quality control after removing low quality reads
  cat(paste0("fastqc -o qa/umitools",
             " umitools/umi_",as.vector(samples$SampleName)[i],"_1.fastq.gz;\n"))
  
  
  #remove low quality reads
  cat(paste0("TrimmomaticSE -threads ",threads,
             " -phred33",
             " umitools/umi_",as.vector(samples$SampleName)[i],"_1.fastq.gz",
             " trimmomatic/TRIM_",as.vector(samples$SampleName)[i],"_1.fastq",
             " HEADCROP:",HEADCROP," ILLUMINACLIP:",ILLUMINACLIP,
             " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36",
             " 2> log/",as.vector(samples$SampleName)[i],"_trim.log;\n"))
  
  
  #quality control after removing low quality reads
  cat(paste0("fastqc -o qa/trimmomatic",
              " trimmomatic/TRIM_",as.vector(samples$SampleName)[i],"_1.fastq;\n"))
  
  #remove rRNA reads
  cat(paste0("time sortmerna ",ref,
              " --reads trimmomatic/TRIM_",as.vector(samples$SampleName)[i],"_1.fastq",
              " --aligned ",as.vector(samples$SampleName)[i],"_aligned",
              " --other ",as.vector(samples$SampleName)[i],"_sortmerna",
              " --log -a ",threads," -v --fastx;\n"))
  
  
  #cmd5_1=paste0("sed '/^$/d' Hif1_D_A_sortmerna.fa > repaired.fa")
  cat(paste0("mv ",as.vector(samples$SampleName)[i],"_sortmerna.fastq",
              " sortmerna/sort_",as.vector(samples$SampleName)[i],"_1.fastq;\n"))
  

  cat(paste0("pigz -p 12 trimmomatic/TRIM_",as.vector(samples$SampleName)[i],"_1.fastq;\n",
              "pigz -p 12 ",as.vector(samples$SampleName)[i],"_aligned.fastq;\n",
              "mv ",as.vector(samples$SampleName)[i],"_aligned.fastq.gz sortmerna/;\n",
              "mv ",as.vector(samples$SampleName)[i],"_aligned.log log/;\n"))
  
  #quality control before alignment
  cat(paste0("fastqc -o qa/sortmerna",
              " sortmerna/sort_",as.vector(samples$SampleName)[i],"_1.fastq;\n"))
  
  #alignment
  cat(paste0("mkdir STAR/",as.vector(samples$SampleName)[i],";\n"))
  
  cat(paste0("~/STAR-2.7.0f/bin/Linux_x86_64/STAR",
              " --genomeDir ",genomeDir,
              " --readFilesIn",
              " sortmerna/sort_",as.vector(samples$SampleName)[i],"_1.fastq",
              " --runThreadN ",threads,
              " --outFileNamePrefix STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],
              " --outFilterMultimapNmax 1",
              " --outSAMtype BAM SortedByCoordinate;\n"))
  

  cat(paste0("pigz -p 12 sortmerna/sort_",as.vector(samples$SampleName)[i],"_1.fastq;\n"))
  
  
  cat(paste0("cp STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],"Log.final.out log/",as.vector(samples$SampleName)[i],"_STAR.log;\n"))
  
  #reads counting
  cat(paste0("htseq-count -f bam -i ",ID,
               " -r pos -m union -s yes",
               " STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],"Aligned.sortedByCoord.out.bam",
               " ",GTF," > count/",as.vector(samples$SampleName)[i],"_union_name.count;\n"))
  
  #contamination test
  cat(paste0("~/fastq_screen_v0.11.1/fastq_screen --aligner bwa --threads ",threads," --conf /mnt/f/Genomes/fastq_screen/fastq_screen.conf ",
               "sortmerna/sort_",as.vector(samples$SampleName)[i],"_1.fastq.gz;\n"))
  cat(paste0("mv sort_",as.vector(samples$SampleName)[i],"_1_screen* qa/fastq_screen;\n"))
  
  
  
  #############UMI TOOLS######
  
 
  #index BAM file and dedup UMI
  
  cat(paste0("samtools index",
               " STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],"Aligned.sortedByCoord.out.bam;\n"))
  
  cat(paste0("umi_tools dedup -I",
               " STAR/",as.vector(samples$SampleName)[i],"/",as.vector(samples$SampleName)[i],"Aligned.sortedByCoord.out.bam",
               " --output-stats=log/dedup/",as.vector(samples$SampleName)[i],
               " -S umitools/",as.vector(samples$SampleName)[i],"-deduplicated.bam",
               " > log/",as.vector(samples$SampleName)[i],"_dedup.log;\n"))
  
  
  #reads counting
  cat(paste0("htseq-count -f bam -i ",ID,
               " -r pos -m union -s yes",
               " umitools/",as.vector(samples$SampleName)[i],"-deduplicated.bam",
               " ",GTF," > count_umi/",as.vector(samples$SampleName)[i],"_union_name.count;\n"))
}

