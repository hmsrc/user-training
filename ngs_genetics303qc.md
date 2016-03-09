# HMS Research Commputing: Intro to Next-Gen Sequencing Technologies Part II

This is the practicum portion of "Intro to Next-Gen Sequencing Technologies," Spring 2016 for HMS Genetics 303qc.  It provices an Orchestra HPC environment-oriented workflow for RNA-seq analysis.

### Logging into Orchestra

Orchestra is the Harvard Medical School High-Performance Compute environment.  With over 7000 cores, 40TB of ram and 25PB attached netework storage, Orchestra is designed to meet the demands of next-generation sequence analysis.  Orchestra allows users to leverage the power of compute across multiple cores in a highly-scalable manner.

Orchestra user accounts can be created for anyone with an eCommons ID.  To have an account created, visit
http://rc.hms.harvard.edu/#orchestra
and fill out the account request form.

To log into Orchestra, 

   * Mac: from the terminal, type
   ```sh
   ssh -l userid orchestra.med.harvard.edu
   ```
   * Windows: from PuTTY, type in  "Host Name"
   ```sh
   orchestra.med.harvard.edu
   ```
   and use your eCommons ID and password
   * Linux: from the terminal, type
   ```sh
   ssh -l userid orchestra.med.harvard.edu
   ```

### Welcome!

Now on the shell login servers (loge/mezzanine), please don't do computationally-heavy work here!
Instead, launch an interactive session to work on a compute node (clarinets, bassoons, fifes, piccolos, etc)

```sh
mfk8@loge:~$ bsub -Is -q interactive bash
mfk8@clarinet002:~$
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

### Orchestra Jobs

Users do their computations on the cluster by submitting jobs (bsubs) to LSF, the cluster scheduler.  LSF identifies the resources needed and queues up jobs, giving priority with higher fairshare (they have submitted less jobs), and to queues with shorter resource requirements.  

Jobs are submitted with a "bsub", and always require
* -q a queue #see table below
* -W a wall time #jobs are killed once they exceed this limit

| Queue Name | Max Cores | Max Runtime |
| ---------- | --------- | ----------- |
| Interactive| 12 | 12hr |
| priority | 20 | 1 month |
| mcore | 20 | 1 montn |
| mpi | 512 | 1 month |
| short | 12 | 12 hr |
| mini | 12 | 10min |
| long | 12 | 1 month |

additional options include
* -n number of course requested #most compute nodes have 12 cores, a small subset has 20
* -R "select[mem=16000]" memory request #most machines have 90GB free, a small subset has 120GB
* -e errorfile
* -o outfile
* -N job completion notification

### Software on Orchestra

Orchestra uses environment modules to manage software and add the appropriate files to a user's PATH.  Environment modules source the most current directories for programs, and solve co-dependencies by loading dependent software.  

```sh
mfk8@clarinet002:~$ module avail                    #see all software available
mfk8@clarinet002:~$ module avail seq                #see only seq software
mfk8@clarinet002:~$ module load seq/fastqc/0.11.3   #load software into your environment
mfk8@clarinet002:~$ module list                     #list all currently loaded modules
mfk8@clarinet002:~$ module unload seq/fastqc/0.11.3 #unload software
mfk8@clarinet002:~$ module list                     #list all currently loaded modules
mfk8@clarinet002:~$ module purge                    #unlod all modules
```

### Monitoring Jobs on Orchestra
```sh

mfk8@clarinet002:~$ bjobs           #lists all jobs running/pending
mfk8@clarinet002:~$ bjobs -l jobid  #gives command used for job
mfk8@clarinet002:~$ bkill jobid     #kills job
mfk8@clarinet002:~$ bkill 0         #kills all jobs
```

### Getting Data To/From Orchestra

You can use an sFTP client to download files to your laptop/desktop. RC recommends "FileZilla," which works on all platforms. Login is the same, orchestra.med.harvard.edu , and files can be dragged and dropped to-from Orchestra.  Simple directory manipulations can also be performed via the GUI, but a "refresh" is required to see the effects.

https://filezilla-project.org/

### Downloading data from GEO

The GeneOmnibus Respository contains a wealth of 'seq experiments, stored as SRA archives.  With the tool "sratoolkit", these can be downloaded and extracted into the native .fastq format, using the command "fastq-dump".  A sample command for downloading an SRA, as a job, looks like this:

```sh
$ module load seq/sratoolkit/2.5.2
$ bsub -q short -W 1:00 "fastq-dump SRRXXXXX"
```

# RNA-seq processing exercise

We will be working with a small toy Drosopila dataset from GEO to familiarize you with an RNA-seq processing workflow.  We will run quality control analysis via FastQC to identify any issues with the runs.  Then we will align these files to the dm3 genome with a popular aligner, TopHat. We will practice manipulating the files using Samtools.  We will download our aligned files to a personal machine, and visuzlize using IGV.  We will then collapse the alignments into counts files using HTSeq, to prepare them for downstream differential expression analysis using tools available in R.

### Class datafiles

Create a directory called "ngsclass", change to it, and copy the class data files from /groups/rc-training/ngsclass to this directory

```sh
mfk8@clarinet002:~$ mkdir ngsclass
mfk8@clarinet002:~$ cd ngsclass
mfk8@clarinet002:~/ngsclass$ cp /groups/rc-training/ngsclass/* .
```

What is your data?

```sh 

mfk8@clarinet002:~/ngsclass$ ls -lh
-rw-rw-r-- 1 kmh40 rccg 101288017 Feb 26 13:01 g1_s1_1.fastq
-rw-rw-r-- 1 kmh40 rccg 101288017 Feb 26 13:01 g1_s1_2.fastq
-rw-rw-r-- 1 kmh40 rccg 102040204 Feb 26 13:01 g1_s2_1.fastq
-rw-rw-r-- 1 kmh40 rccg 101959796 Feb 26 13:01 g1_s2_2.fastq
-rw-rw-r-- 1 kmh40 rccg 101540796 Feb 26 13:01 g2_s1_1.fastq
-rw-rw-r-- 1 kmh40 rccg 101459204 Feb 26 13:01 g2_s1_2.fastq
-rw-rw-r-- 1 kmh40 rccg 101288017 Feb 26 13:01 g2_s2_1.fastq
-rw-rw-r-- 1 kmh40 rccg 101197361 Feb 26 13:01 g2_s2_2.fastq

```

This is a Drosophila fastq dataset, with 2 groups (g1/g2) with 2 samples each (s1/21), paired end (_1, _2).  But really, what is your data?

```sh
$ head g1_s1_1.fastq
@SRR2107137.1 HISEQ:202:C5R67ACXX:2:1101:1246:1937 length=202
AAAATTACTTTTTATTAACCGTGTTTTCTTGGCCGTCAATTACGTTGCTTTTCTTTTCGCTGCGTTTCCGGAAGAGGATATTGAGTAGCACCGTTTCGGGCCCTTGGTCTTGGTCATTGTGCTGTTTAGCCCGAAACGGTGCTACTCAATATCCTCTTCCGGAAACGCAGCGAAAAGAAAAGCAACGTAATTGACGGCCAAG
...
```
The SRA identity is found in the @ line.  You can search GEO for this SRR number to get more information on the experiment.


### QC Check: FastQC

A quick QC check will identify problems with the quality of the reads, identify adapter/barcode sequence, kmers, and more.  FastQC is the standard for performing efficient QC checks.  It creates an html report for each file.  These html reports are best downloaded and viewed on  a personal computer.  

We will submit jobs to Orchestra to perform FastQC.  This program will take less than 5 minutes to complete, and only requires 1 core, so the jobs belong in the "short" queue.  However, since we only have a few jobs to run, and 2 can run at a time in the "priority" queue with little to no pend time, we will run our jobs there.  With a simple "for" loop, we can submit all of the jobs at once, by looping through all of the files ending in ".fastq" in the directory.

```sh
$ module load seq/fastqc/0.11.3 
$ for i in *.fastq; do bsub -q priority -W 5 fastqc $i; done
```

We will download these files to our personal computers to view them.  In FileZilla, navigate to your "seqclass" folder, and drag and drop the files labeled below to a location on your computer.
```sh
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

These files are 89-91bp long, with Sanger 1.5 PHRED encoding.  They are of acceptable quality, with good GC content and few repeated sequences.  These files are not good candidates for trimming, so we can proceed with the alignment.

### TopHat2 Alignment

We will first use TopHat2 to align these files to a reference genome,  and annotate them.  TopHat2 relies on its partner program Bowtie1/2 to do the alignment, while TopHat2 maps splice junctions, transcription splice sites, and novel isoforms.  

We will first load the TopHat2/Bowtie2 modules, along with Samtools.  Note all of the other software co-loaded by Orchestra, to solve dependency issues.

```sh
$ module load seq/tophat/2.1.0 seq/bowtie/2.1.0 seq/samtools/1.2
$ module list
Currently Loaded Modulefiles:
  1) dev/java/jdk1.8          4) atlas/3.10.2             7) utils/hdf5/1.8.15       10) seq/samtools/1.2
  2) seq/fastqc/0.11.3        5) dev/compiler/gcc-4.8.5   8) dev/python/2.7.6
  3) seq/bowtie/2.1.0         6) utils/zlib/1.2.8         9) seq/tophat/2.1.0


```

Bowtie relies on index files to speedily align to a reference genome.  These are Bowtie-parsed versions of the reference genome in a format that Bowtie can read.  Bowtie1 index files end in .ebwt, and Bowtie2 index files end in .bt2 .  We are using Bowtie2.  The Illumina igenome project created Bowtie index files for common genomes to standard reference genomes.  In Orchestra, these are located in /groups/shared_databases/igenomes/organism .  For Bowtie2 indexes for hg19, hg18, mm10, mm9, and dm3, we have created a softlink that points directily to these igenome Bowtie2 indexes that are located in /groups/shared_databases/igenome/organism/UCSC/version/Sequence/Bowtie2Index, and can just be referenced as such.  

One of the key differences between UCSC and NCBI notation is how chromosomes are called.  In UCSC, the chromosome is called by "chr1", in Ensembl, it is just "1".  In order to use GTF annotation files on UCSC-aligned .bam files, the "chr" must first be stripped from the GTF file.

For this TopHat2 alignment, we are considering the sequencing library prep to be unstranded.

We will be utilizing multithreading to distribute the compute job over multiple cores (the -p option).  The majority of Orchestra machines have up to 12 cores available per node; a small subset have 20.  The more cores that are requested, the longer a job takes to dispatch, as resources are collected for the job.  

The TopHat command format is: tophat -p processors -o outputdirectory path/to/genomeIndexFiles read1.fastq (read2.fastq)

```sh
$ mkdir tophat_logs
$ mkdir tophat_g1_s1
$ bsub -q mcore -W 1:00 -n 4 -o ./tophat_logs/%J.out -e ./tophat_logs/%J.err -N "tophat -p 4 -o ./tophat_g1_s1 dm3 g1_s1_1.fastq g1_s1_2.fastq"
$ mkdir tophat_g1_s2
$ bsub -q mcore -W 1:00 -n 4 -o ./tophat_logs/%J.out -e ./tophat_logs/%J.err -N "tophat -p 4 -o ./tophat_g1_s2 dm3 g1_s2_1.fastq g1_s2_2.fastq"
$ mkdir tophat_g2_s1
$ bsub -q mcore -W 1:00 -n 4 -o ./tophat_logs/%J.out -e ./tophat_logs/%J.err -N "tophat -p 4 -o ./tophat_g2_s1 dm3 g2_s1_1.fastq g2_s1_2.fastq"
$ mkdir tophat_g2_s2
$ bsub -q mcore -W 1:00 -n 4 -o ./tophat_logs/%J.out -e ./tophat_logs/%J.err -N "tophat -p 4 -o ./tophat_g2_s2 dm3 g2_s2_1.fastq g2_s2_2.fastq"
```

Where are your aligned files?  For each analysis, you created a folder called tophat_gX_sX, and the aligned file is called "accepted_hits.bam".  The summary of the alignment is "align_summary.txt:, and the insertions, deletions, and splice junctions are bed coordinates files named accordingly.  Information on job runtime/success/errors is found in the folder "tophat_logs", by jobid.  Look at these files by "less tophat_logs/JOBID.err" to view.

Let's rename the .bam files and create .bai (.bam index files) for each TopHat hit.  These will be necessary for visualization.

```sh
# from root data directory
$ cp  tophat_g1_s1/accepted_hits.bam tophat_g1_s1/tophat_g1_s1.bam
$ bsub -q short -W 5 "samtools index tophat_g1_s1/tophat_g1_s1.bam"
$ cp tophat_g1_s2/accepted_hits.bam tophat_g1_s2/tophat_g1_s2.bam
$ bsub -q short -W 5 "samtools index tophat_g1_s2/tophat_g1_s2.bam"
$ cp  tophat_g2_s1/accepted_hits.bam tophat_g2_s1/tophat_g2_s1.bam
$ bsub -q short -W 5 "samtools index tophat_g2_s1/tophat_g2_s1.bam"
$ cp tophat_g2_s2/accepted_hits.bam tophat_g2_s2/tophat_g2_s2.bam
$ bsub -q short -W 5 "samtools index tophat_g2_s2/tophat_g2_s2.bam"

```

How did the aligment perform?  TopHat creates summary statistics for each run, called "align_summary.txt".

```sh
# Here is the reported statistics for TopHat, g1_s1

$ less tophat_g1_s1/align_summary.txt
Left reads:
          Input     :    500000
           Mapped   :    500000 (100.0% of input)
            of these:      5879 ( 1.2%) have multiple alignments (3374 have >20)
Right reads:
          Input     :    500000
           Mapped   :    500000 (100.0% of input)
            of these:      5879 ( 1.2%) have multiple alignments (3374 have >20)
100.0% overall read mapping rate.

Aligned pairs:    500000
     of these:      5879 ( 1.2%) have multiple alignments
                  500000 (100.0%) are discordant alignments
 0.0% concordant pair alignment rate.
```

### Read Visualization with IGV

BAM files and BED files can be visualized with the Java-based tool IGV.  IGV can be launched with larger instances of memory by modifying the .bat file; for this exercise, the default is fine. 

First download the BAM files with their corresponding BAI folders to your computer, using FileZilla.  IGV needs both the .bam and its associated .bai to load properly.

Then, launch IGV, and load the reference genome track (Genomes->Load Genome from Server->D. melanogaster (dm3)).  Now, you can load your BAM files (File->Load from File)



### Read Counting with HTSeq

The number of reads assigned to a gene feature can be turned into a counts file with htseq-count.  This Python tool relies on the GTF coordinates and annotation to assign a read to a gene.  Files must first be name sorted with Samtools to create a "*.sort.bam file" and can then be counted. We are assuming the experiments are not stranded.

```sh
# Sort the TopHat experiment on name (-n)
$ bsub -q short -W 5 "samtools sort -n tophat_g1_s1/tophat_g1_s1.bam tophat_g1_s1/tophat_g1_s1.sort"
$ bsub -q short -W 5 "samtools sort -n tophat_g1_s2/tophat_g1_s2.bam tophat_g1_s2/tophat_g1_s2.sort"
$ bsub -q short -W 5 "samtools sort -n tophat_g2_s1/tophat_g2_s1.bam tophat_g2_s1/tophat_g2_s1.sort"
$ bsub -q short -W 5 "samtools sort -n tophat_g2_s2/tophat_g2_s2.bam tophat_g2_s2/tophat_g2_s2.sort"
```

```sh
# Load htseq module
$ module load seq/htseq/0.6.1

# command reads: htseq-count options sortedBamFile.bam /path/to/GTF.gtf > output.txt

$ bsub -q short -W 30 "htseq-count  --order=name --stranded=no --format=bam tophat_g1_s1/tophat_g1_s1.sort.bam /groups/shared_databases/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf > tophat_g1_s1/tophat_g1_s1.counts.txt"
$ bsub -q short -W 30 "htseq-count  --order=name --stranded=no --format=bam tophat_g1_s2/tophat_g1_s2.sort.bam /groups/shared_databases/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf > tophat_g1_s2/tophat_g1_s2.counts.txt"
$ bsub -q short -W 30 "htseq-count  --order=name --stranded=no --format=bam tophat_g2_s1/tophat_g2_s1.sort.bam /groups/shared_databases/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf > tophat_g2_s1/tophat_g2_s1.counts.txt"
$ bsub -q short -W 30 "htseq-count  --order=name --stranded=no --format=bam tophat_g2_s2/tophat_g2_s2.sort.bam /groups/shared_databases/igenome/Drosophila_melanogaster/UCSC/dm3_2L_2R/Annotation/Genes/genes.gtf > tophat_g2_s2/tophat_g2_s2.counts.txt"
```

We can now verify that the counting worked by looking at the "head" and "tail" of our counts.txt files .  Note that the uncounted read statistics comprise the last 5 lines of this file.  These may need to be trimmed out in later downstream analyses.

```sh
$ head tophat_g1_s1/tophat_g1_s1.counts.txt
128up   151
14-3-3zeta      486
18w     47
5-HT1A  1
5-HT1B  0
A16     117
ACXA    0
ACXB    0
ACXC    0
ACXE    0

$ tail tophat_g1_s1/tophat_g1_s1.counts.txt

yuri    9
zetaTry 13
zf30C   207
zip     415
zuc     31
__no_feature    5734
__ambiguous     8461
__too_low_aQual 0
__not_aligned   0
__alignment_not_unique  81056

```
### Next Steps

These read counts files can be aggregated in bash or R and run through popular R-based differential expression algorithms like edgeR and DESeq2, and used for sophisticated plotting.  The combined counts matrices can even be imported into programs like MS Excel for filtering, low-level statistics, and graphing.  It is of note that these counts files have not been normalized in any way; R-based Differential Expression algorithms control for effects based on library size, not gene length.  Gene length can be accounted for in functional enrichment analysis, using tools like the R-based GOSeq.  