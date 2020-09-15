# HMS Research Computing: Intro to Next-Gen Sequencing Technologies Part II

This is the practicum portion of "Intro to Next-Gen Sequencing Technologies," Spring 2018 for HMS Genetics 303qc.  It provides an O2 HPC environment-oriented workflow for RNA-seq analysis.

### Logging into O2

O2 is the Harvard Medical School High-Performance Compute environment.  With over 8000 cores, 45 TB of RAM and 32 PB attached network storage, O2 is designed to meet the demands of next-generation sequence analysis.  O2 allows users to leverage the power of compute across multiple cores in a highly scalable manner.

O2 user accounts can be created for anyone with an eCommons ID.  To have an account created, visit
http://rc.hms.harvard.edu/#cluster
and fill out the account request form.

To log into O2, 

   * Mac: from the terminal, type
   ```
   ssh yourEcommons@o2.hms.harvard.edu
   ```
   * Windows: from MobaXterm, type
   ```
   ssh yourEcommons@o2.hms.harvard.edu
   ```
   * Linux: from the terminal, type
   ```
   ssh yourEcommons@o2.hms.harvard.edu
   ```

### Welcome!

Now you're on the shell login servers (login01-05), please don't do computationally-heavy work here!
Instead, launch an interactive session to work on a compute node (compute-a, compute-e etc)

```sh
mfk8@login01:~$ srun --pty -p interactive -t 0-12:00 --mem 8G bash
mfk8@compute-a:~$
```

### Class datafiles

Create a directory in /n/scratch2 (10TB limit) called "ngsclass", change to it, and copy the class data files from /n/groups/rc-training/ngsclass to this directory

```sh
mfk8@compute-a:~$ mkdir -p /n/scratch2/$USER/ngsclass               #create folder "ngsclass" in /n/scratch2 under your user name
mfk8@compute-a:~$ cd /n/scratch2/$USER/ngsclass                     #change directory to ngsclass
mfk8@compute-a:/n/scratch2/mfk8/ngsclass$ cp -r /n/scratch2/rc-training/ngsclass/* . #download all contents to current location (.)
```

### Linux 101

Cheatsheet of commands
```sh
cd = Change Directory (".." is up one level, "." is present level)
cp = CoPy (from, to)
ls = LiSt directory contents
ls -l = LiSt Long form: get file details
less = read flat text file contents, scroll up/down, "q" to quit
mkdir = MaKe DIRectory
head = Display the first 10 lines
tail = Display the last 10 lines
```

### O2 Jobs

Users do their computations on the cluster by submitting jobs (sbatchs) to SLURM, the cluster scheduler.  SLURM identifies the resources needed and queues up jobs, giving priority to users with higher fairshare (they have submitted less jobs), to partitions with shorter time limits, and to jobs with smaller resource requirements.  

Jobs are submitted with a "sbatch", and always require
* -p a partition #see table below
* -t a wall time #jobs are killed once they exceed this limit

| Queue Name | Max Cores | Max Runtime | Best for |
| ---------- | --------- | ----------- | -------- |
| interactive| 20 | 12hr | programming, interactive use |
| priority | 20 | 1 month | urgent jobs (2 max) |
| short | 20 | 12 hr | short single/multicore jobs |
| medium | 20 | 5 days | medium single/multicore jobs |
| long | 20 | 1 month | long single/multicore jobs |
| mpi | 640 | 5 days | large jobs with special parallel code |
| gpu | - | 72 gpu hours | GPU jobs |
| transfer | 4 | 5 days | Moving files to O2 from research.files |
| highmem | 1TB | 5 days | Large memory jobs |

additional options include
* -c number of course requested #max 20
* --mem  memory request #most machines have 250GB available
* -e %j.err errorfile with jobid inserted
* -o %j.out outfile with  jobid inserted

### Software on O2

O2 uses LMOD modules to manage software and add the appropriate files to a user's $PATH.  Environment modules source the most current directories for programs, and solve co-dependencies by loading dependent software.  

```sh
mfk8@compute-a:~$ module spider                         #see all software available
mfk8@compute-a:~$ module load gcc/6.2.0 star/2.5.2b     #load software into your environment
mfk8@compute-a:~$ module list                           #list all currently loaded modules
mfk8@compute-a:~$ module unload star/2.5.2b             #unload software
mfk8@compute-a:~$ module purge                          #unload all modules
```

### Monitoring Jobs on Orchestra
```sh

mfk8@compute-a:~$ squeue -u $USER                                         #lists all jobs I have running/pending
mfk8@compute-a:~$ sacct -u $USER --format=JobID,JobName,MaxRSS,Elapsed    #accounting details
mfk8@compute-a:~$ scancel jobid                                           #kills job
```

### Getting Data To/From Orchestra

You can use an sFTP client to download files to your laptop/desktop. RC recommends "FileZilla," which works on all platforms. Login is to `transfer.rc.hms.harvard.edu`, port `22` and files can be dragged and dropped to-from O2.  Simple directory manipulations can also be performed via the GUI, but a "refresh" is required to see the effects.

https://filezilla-project.org/

### Downloading data from GEO

The GeneOmnibus Respository contains a wealth of 'seq experiments, stored as SRA archives.  With the tool "sratoolkit", these can be downloaded and extracted into the native .fastq format, using the command "fastq-dump".  A sample command for downloading an SRA, as a job, looks like this:

```sh
$ less scripts/getSRA.run
```

```sh
#!/bin/bash                     
#SBATCH -p short                #short partition
#SBATCH -t 0-1:00               #1 hour time limit
#SBATCH -e %j.err               #error log
#SBATCH -o %j.out               #out log
module load sratoolkit/2.8.1    #load module
fastq-dump --split-files $1     #run fastq-dump (split paired end fastq) on command line input file name
```
executed as

```sh
$ sbatch scripts/getSRA.run SRAnumberhere
```

# RNA-seq processing exercise

We will be working with a small, unpublished toy mouse dataset from GEO to familiarize you with an RNA-seq processing workflow.  We will run quality control analysis via FastQC to identify any issues with the runs, and collate the reports with MultiQC.  Then we will align these files to the NCBI GRCm38 genome with a popular aligner, STAR. We will practice manipulating the files using Samtools.  For visualization of the reads using IGV, we will download our aligned files to a personal machine. Last, we will use the counts files to perform differential expression analysis using the edgeR packgage in R.

### Inspecting data

What is your data?

```sh 

mfk8@compute-a:/n/scratch2/mfk8/ngsclass$ ls -lh
total 66G
-rw-rw-r-- 1 kmh40 rccg 6.9G Mar 13 15:06 g1_s1_1.fastq
-rw-rw-r-- 1 kmh40 rccg 6.9G Mar 13 15:08 g1_s1_2.fastq
-rw-rw-r-- 1 kmh40 rccg 5.9G Mar 13 15:10 g1_s2_1.fastq
-rw-rw-r-- 1 kmh40 rccg 5.9G Mar 13 15:11 g1_s2_2.fastq
-rw-rw-r-- 1 kmh40 rccg 7.6G Mar 13 15:13 g2_s1_1.fastq
-rw-rw-r-- 1 kmh40 rccg 7.6G Mar 13 15:14 g2_s1_2.fastq
-rw-rw-r-- 1 kmh40 rccg 7.4G Mar 13 15:16 g2_s2_1.fastq
-rw-rw-r-- 1 kmh40 rccg 7.4G Mar 13 15:16 g2_s2_2.fastq
drwxrwsr-x 2 kmh40 rccg  516 Mar 13 17:42 GRCm38_star_125
drwxrwsr-x 2 kmh40 rccg  150 Mar 13 17:46 scripts
```

This is a Mouse fastq dataset, with 2 groups (g1/g2) with 2 samples each (s1/s2), paired end (_1, _2).  But really, what is your data?

```sh
$ head g1_s1_1.fastq
@SRR6725731.1 HISEQ:473:CB152ANXX:8:1101:2974:1993 length=126
NTCTAGAGCTAATACATGCCGACGGGCGCTGACCCCCCTTCCCGGGGGGGGATGCGTGCATTTATCAGCCAGATCGGAAGGGCACACGTCTGAACTCCAGTCACCGTCTAACCACTACAGATCTCG
+SRR6725731.1 HISEQ:473:CB152ANXX:8:1101:2974:1993 length=126
#<=@BFGGGEGGFGGGGGGGGGGGFGGGGGGGGGGGGGEBGGGGGGGGGGGGGGGGGGGGGD=EGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG9/DGG
...
```

The SRA identity is found in the @ line.  You can search GEO for this SRR number to get more information on the experiment.


### QC Check: FastQC

A quick QC check will identify problems with the quality of the reads, identify adapter/barcode sequence, kmers, and more.  FastQC is the standard for performing efficient QC checks.  It creates an html report for each file.  We will aggregate the reports together using MultiQC.  These html reports are best downloaded and viewed on a personal computer.  

We will submit jobs to O2 to perform FastQC.  This program will take less than 5 minutes to complete, and only requires 1 core with 8GB of memory, so the jobs belong in the "short" queue. With a simple "for" loop, we can submit all of the jobs at once, by looping through all of the files ending in ".fastq" in the directory.

The sbatch scripts for these lessons are contained in your "scripts" folder.  To interface with the scheduler, we wrap up our commands into a script, and submit that script to the scheduler.  These scripts can take additional command line arguments (like file names, or options), and can be passed after the script name in numerical order as $1, $2, etc.

The contents of fastQC.run can be viewed with

```sh
$ less scripts/fastQC.run
```

```
#!/bin/bash                 
#SBATCH -p short            #partition
#SBATCH -t 0-00:05          #time limit, 5min
#SBATCH --mem 1G            #memory requested
module load fastqc/0.11.3   #software to use
fastqc $1                   #command
```
We sill execute this with a for loop saying for every file ending in .fastq, run fastqQC on it.

```sh
$ for file in *.fastq; do sbatch scripts/fastQC.run $file; done
```
We can check on these jobs status by running

```sh
$ squeue -u $USER
```

Once these jobs are all finished running, we can aggregate the reports together with MultiQC.  This takes our 8 `fastq.html` reports, and we will create a master report `multiqc_report.html`.  Since we're in an interactive session (terminal says compute-a or compute-e), we can load the MultiQC module directly and run, it takes less than a minute

```sh
$ module load gcc/6.2.0 python/2.7.12 multiqc/1.3
$ multiqc . #runs multiqc on all fastqc files in the directory
```
We will download these files to our personal computers to view them.  In FileZilla, navigate to your "ngsclass" folder, and drag and drop the files labeled below to a location on your computer.
```sh
multiqc_report.html
g1_s1_1.fastqc.html
g1_s1_2.fastqc.html
g1_s2_1.fastqc.html
g1_s2_2.fastqc.html
g2_s1_1.fastqc.html
g2_s1_2.fastqc.html
g2_s2_1.fastqc.html
g2_s2_2.fastqc.html
```

To view, open your file exploring manager and double-click on the .html report.  The report will open in your default browser.

These files are 126bp long, with Sanger 1.9 PHRED encoding.  They are of acceptable quality, with good GC content and few repeated sequences.  There is a lot of adapter content in the sequences.  These files may benefit from adapter trimming, but our aligner will soft clip (ignore) the 5' and 3' ends to align, and we will have an acceptable map rate.

### STAR Alignment

We will first use STAR to align these files to a reference genome,  and count the reads assigned to each gene.  STAR will create a `.bam` file, which is how the reads mapped, including mapping quality, and a counts file, which we will use to run the differential expression analysis.  

STAR relies on index files to speedily align to a reference genome.  These are STAR-parsed versions of the reference genome in a format that STAR can read.  Currently in O2, there are some STAR indices available, but we recommend creating your own using STAR genome-generate on the FASTAs (genomic sequence) and GTF file (annotation file) from the organism and build (UCSC, NCBI) that you choose. Here, the NCBI GRCm38, version 91 STAR indices with a splice junction overhang 125 (read length-1) were created and will be referenced out the "GRCm38_star_125" folder you downloaded. 

One of the key differences between UCSC and NCBI notation is how chromosomes are called.  In UCSC, the chromosome is called by "chr1", in NCBI/Ensembl, it is just "1".  In order to use GTF annotation from one format on the other, the "chr" must be added or deleted.

For this STAR alignment, we are considering the sequencing library prep to be unstranded.

We will be utilizing multithreading to distribute the compute job over multiple cores (the --runThreadN option).  The majority of O2 machines have up to 32 cores available per node; 20 are permitted to be used for each multicore job.  The larger number of cores that are requested, the longer a job takes to dispatch, as resources are collected for the job. Performance does not scale linearly; the more cores requested does not mean that much speedup, many NGS algorithms top out at usage of 6-8 cores, and just require sufficient memory.  STAR typically just requires 2 cores and a fair bit (48G) of memory.

The STAR  command format is: 
```sh
STAR 
--runThreadN #numer of cores 
--genomeDir #where the genome indices are located
--sjdbGTFfile #where the GTF file is located
--readFilesIn read1.fastq read2.fastq #for non-zipped files
--outFileNamePrefix "${1%.*}"_star #take Read 1, strip the suffix, and name all output files with this convention + "star" 
--outSAMtype BAM SortedByCoordinate #write an aligned BAM file and sort in chromosomal order
--outSAMunmapped Within #retain unmapped reads in the BAM file
--outSAMattributes NH HI NM MD AS #extra BAM file attributes
--outReadsUnmapped Fastx #write unmapped reads to a FASTX, can be used for blastn/downstream triage
--quantMode GeneCounts #count reads assigned to each gene
```

We will be mirroring the resource request options in our sbatch (-runThreadN 2 cores requested in STAR, -c 2 cores requested in O2 bsub).  As a good practice, we are creating error and output files, where %j inserts the jobid that O2 allocates.  These error and output files (logs) will be collected in your working directory.

The STAR scripts's contents are as follows:

```sh
$ less scripts/star.run
```

```sh
#!/bin/bash
#SBATCH -p short               #partition
#SBATCH -t 1:00:00             #wall time
#SBATCH -c 2                   #cores requested
#SBATCH --mem 48G              #memory requested
#SBATCH -o star.%j.out         #job out logs
#SBATCH -e star.%j.err         #job error logs

module load gcc/6.2.0 star/2.5.2b
STAR --runThreadN 2 --genomeDir  GRCm38_star_125 --sjdbGTFfile GRCm38_star_125/Mus_musculus.GRCm38.91.chr.gtf --readFilesIn $1 $2 --outFileNamePrefix "${1%.*}"_star --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes NH HI NM MD AS --outReadsUnmapped Fastx --quantMode GeneCounts
```

We will execute this script on each pair of reads (_1 and _2) by passing them as command line arguments ($1 and $2) to the sbatch script.

```sh
$ sbatch scripts/star.run g1_s1_1.fastq g1_s1_2.fastq #group 1 sample 1 reads 1 and 2
$ sbatch scripts/star.run g1_s2_1.fastq g1_s2_2.fastq #group 1 sample 2 reads 1 and 2
$ sbatch scripts/star.run g2_s1_1.fastq g2_s1_2.fastq #group 2 sample 1 reads 1 and 2
$ sbatch scripts/star.run g2_s2_1.fastq g2_s2_2.fastq #group 2 sample 2 reads 1 and 2
```

We can monitor these jobs again with our "squeue" command
```sh
$ squeue -u $USER
```

Where are your aligned files?  For each pair of reads, STAR creates a .bam file (aligned reads) and ReadsPerGene.out.tab (counts file).  
How did the alignment perform?  The summary of the alignment is "starLog.final.out".

Let's take a look at the summary statistics for group 1, sample 1:

```sh
$ less g1_s1_1_starLog.final.out
```

```sh
                                 Started job on |       Mar 12 10:30:20
                             Started mapping on |       Mar 12 10:33:36
                                    Finished on |       Mar 12 10:56:01
       Mapping speed, Million of reads per hour |       50.21

                          Number of input reads |       18760811
                      Average input read length |       252
                                    UNIQUE READS:
                   Uniquely mapped reads number |       14108887
                        Uniquely mapped reads % |       75.20%
                          Average mapped length |       244.30
                       Number of splices: Total |       6707959
            Number of splices: Annotated (sjdb) |       6498280
                       Number of splices: GT/AG |       6517075
                       Number of splices: GC/AG |       63195
                       Number of splices: AT/AC |       7736
               Number of splices: Non-canonical |       119953
                      Mismatch rate per base, % |       0.52%
                         Deletion rate per base |       0.11%
                        Deletion average length |       4.09
                        Insertion rate per base |       0.03%
                       Insertion average length |       1.71
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       1056635
             % of reads mapped to multiple loci |       5.63%
        Number of reads mapped to too many loci |       74588
             % of reads mapped to too many loci |       0.40%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       18.16%
                     % of reads unmapped: other |       0.61%
                                  CHIMERIC READS:
                       Number of chimeric reads |       0
                            % of chimeric reads |       0.00%

```

You can also rerun MultiQC on the directory, and it will add the STAR alignment statistics to the report.

### Read Visualization Preparation

If we'd like to visualize our reads, we need to create a BAM index file (`.bai`).  Samtools is a multi-purpose tool for manipulating SAM/BAM files and getting statistics.  Here, we will use samtools to index the aligned BAM files.

```sh
$ less scripts/bamindex.run
```

```sh
#!/bin/bash            
#SBATCH -p priority     #partition
#SBATCH -t 0-00:10      #time, 10 minutes
#SBATCH --mem 8G        #memory in GB
#SBATCH -e bi.%j.err    #error logs
#SBATCH -o bi.%j.out    #out logs
module load gcc/6.2.0 samtools/1.3.1 #load module
samtools index $1                     #index the command line argument
```

We can loop over all of our BAM files, submitting a priority job for each (2 can run at a time).

```sh
$ for bamfile in *.bam; do sbatch scripts/bamindex.run $bamfile; done
```

Once they're finished, we should now see 4 .bai (BAM index) files, one for each BAM file.

```sh
$ ls -lh *.bai
```
### Read Visualization with IGV

BAM files and BED files can be visualized with the Java-based tool IGV.  IGV can be launched with larger instances of memory by modifying the .bat file; for this exercise, the default is fine. 

First download the `.bam` files with their corresponding `.bai` files to your computer, using FileZilla.  IGV needs both the .bam and its associated .bai to load properly.

Then, launch IGV, and load the reference genome track (Genomes->Load Genome from Server->Mouse (mm10)).  Now, you can load your BAM files (File->Load from File)

You'll need to explore navigating IGV by gene names and chromosomal coordinates for the assignment.  Here's a zoomed out view of TP53, referenced as "Trp53" in the navigation bar:

![alt text](https://github.com/hmsrc/user-training/blob/master/trp53.igv.png "IGV Screenshot")


### Counts Files

The counts files are the measure of how many read pairs are assigned to a gene.  There are three types of RNA-seq library prep that can used: unstranded (traditional), stranded, or reverse.  For this library prep, we will consider it to be an "unstranded" prep, since we have no other information.  This is column 2 of each counts file.

**You'll need to use FileZilla to download each `_starReadsPerGene.out.tab` for each read pair in preparation for the next exercises**

### Next Steps

These read counts files can be aggregated read into R and run through popular R-based differential expression algorithms like edgeR and DESeq2, and used for sophisticated plotting.  The combined counts matrices can even be imported into programs like MS Excel for filtering, low-level statistics, and graphing.  It is of note that these counts files have not been normalized in any way; R-based differential expression algorithms control for effects based on library size or median expression, not gene length.  Gene length can be accounted for in functional enrichment analysis, using tools like the R-based GOSeq.  
