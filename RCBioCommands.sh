#!/bin/sh

set -x 

# lool for module
module spider rcbio

# load module
module load rcbio/1.3.3

# insteractive job
srun --pty -p interactive -t 0-2:00  --mem 2G bash

# run bowtie software
bowtie -c 4 hg19 file1_1.fq file1_2.fq

# Bash commands
cd # go to a directory 

ls # list directory and file

cp # copy directory and file

mkdir # create new directory

less # view file content

cat # show file content

# For loop
for i in 1 2 3 ; do mkdir $i ; done

for i in `cat list.txt` ; do cp $i ~/work ; done

# Bash variable
varName=“value”

varName=value

# Refer variable: 
${varName}, $varName, “$varName”, “$varName“ 

# Bash path/file name manipulation

fileWithPath=/home/id/data/tumor.fq

dirname $fileWithPath           # /home/id/data

basename $fileWithPath          # tumor.fq

basename $fileWithPath .fq      # tumor

# Bash string manipulation
someString=abcdefghijk

echo ${someString#abc}            # defghijk

echo ${someString%ijk}             # abcdefgh

prefix=abc; suffix=ijk

echo ${someString#$prefix}      # defghijk
echo ${someString%$suffix}      # abcdefgh

# Bash string concatenate
sample=tumor; rep=2
echo $sample.$ref                  # rumor2
echo $sample_$rep                # does not work          
echo /data/group/$sample/$tumor/$rep

# Catch the output of command
output=`cmd -in $input -out $out`

output=$(cmd -in $input -out $out)

output=$(echo today is $(date))

output=`echo today is \`date\``

# Combine multiple commands
cmd1; cmd2; cmd3;

cmd1 && cmd2 && cmd3

mkdir outputDir && cmdToCreateOutput

ls $input && cmd –in $input || exit 

# Creating a Simple Text File
# “Nano,” “vi”, “emacs” are simple command-line editors available.
# To create a new file, type the editor you want, then the name of the new file.  To edit an existing file, do the same.
nano myfile.txt

# then type this sentence:
		
        This is my new file text.
# (Control-X to save (yes) and exit.)

# check the file list in the directory
ls
# we get	 	
myfile.txt


