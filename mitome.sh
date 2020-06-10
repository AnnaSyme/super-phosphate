#!/usr/bin/env bash

#A script to assemble a mitome

#............................................................................
# How to use

# activate your conda environment with the tools needed
# bash mitome.sh -b baits R1 R2 nano

#............................................................................
# Defaults

set -e #exit if a command exits with non zero status
script=$(basename $0) #script name less file path
baits=""
genome_size=700000
threads=16

#............................................................................
# Functions

function msg {
  echo -e "$*"
}
# $* is args passed in, -e means interpret backslashes

function err {
  echo "Error: $*" 1>&2
}
# redirects stdout to stderr

function banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function msg_banner {
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
  msg "$*"
  printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function usage {
  banner
  msg "$script\n   Assemble a mitome with long and short reads."
  msg "Author:\n  Anna Syme, anna.syme@gmail.com"
  msg "Usage:\n   $script [options] R1 R2 nano"
  msg "Parameters:"
  msg "   Illumina R1 reads, trimmed, filtered      R1.fq.gz"
  msg "   Illumina R2 reads, trimmed, filtered      R2.fq.gz"
  msg "   Nanopore reads, raw                       nano.fq.gz"
  msg "Options:"
  msg "   -h                     Show this help"
  msg "   -t NUM                 Number of threads (default=16)"
  msg "   -g NUM                 genome size in bp (default=700000)"
  msg "   -b FILE                bait sequences file name"
  msg "Example:"
  msg "   $script -t 8 R1.fq.gz R2.fq.gz minion.fq.gz"
  banner
  exit 1
  #exits script
}
#...........................................................................
# Parse the command line options

#loops through, sets variable or runs function
#have to add in a colon after the flag, if it's looking for an arg
#e.g. t: means flag t is set, then look for the arg (eg 32) to set $threads to
# : at the start disables default error handling for this bit
#instead it will label an odd flag as a ?
#the extra : at the end means it will label missing args as :

while getopts ':ht:g:b::' opt ; do
  case $opt in
    h)
      usage
      ;;
    t)
      threads=$OPTARG
      ;;
    g)
      genome_size=$OPTARG
      ;;
    a)
      adapters=$OPTARG
      ;;
    b)
      baits=$OPTARG
      ;;
    \?)
      echo "Invalid option '=$OPTARG'"
      exit 1
      ;;
    :)
      echo "Option '-$OPTARG' requires an argument"
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))
#remove all options that has been parsed by getopts
#so that $1 refers to next arg passed to script
#e.g. a positional arg

if [ $# -ne 3 ]; then
  msg "\n **Please provide three input parameters** \n"
  usage
fi
#if number of pos params is not 3, print usage msg

#...........................................................................
#start
msg "\n"
msg_banner "now running $script"

#...........................................................................
#check inputs

#give variable names to the inputs
R1=$1
R2=$2
nano_raw=$3

msg "This script will use:"
msg "   Illumina reads R1:   $R1"
msg "   Illumina reads R2:   $R2"
msg "   Nanpore reads:       $nano_raw"
msg "   Bait sequences:      $baits"
msg "   Adapter sequences:   $adapters"
msg "   Genome size:         $genome_size"
msg "   Threads:             $threads"

#...........................................................................

conda env export --name bio > bio.yml
#saves conda env with tools and versions

#...........................................................................
msg_banner "Round 1"

#...........................................................................
msg_banner "now extracting mitome nanopore reads from all reads"

nano_mt1=nano_mt1.fq.gz

minimap2 -a -x map-ont -t $threads $baits $nano_raw | 
samtools fastq -0 $nano_mt1 -n -F 4 -

#map the reads to a set of baits (e.g. CDS.fasta from sister taxa) 
#this pipes the output to samtools fastq,
#which extracts the fastq reads from the alignment
#the flag -F 4 means exclude unmapped reads (i.e., non mt reads)

#...........................................................................
msg_banner "now assembling nanopore mt reads"

flye --nano-raw $nano_mt1 --genome-size $genome_size \
--out-dir flye001 --threads $threads

#using uncorrected reads as read correction can lead to errors

assembly1=assembly1.fasta

cp flye001/assembly.fasta $assembly1

#...........................................................................
msg_banner "Round 2"

#...........................................................................
msg_banner "now extracting mitome nanopore reads from all reads"

nano_mt2=nano_mt2.fq.gz

minimap2 -m 5000 -a -x map-ont -t $threads $assembly1 $nano_raw | 
samtools fastq -0 $nano_mt2 -n -F 4 -

#map the reads to a set of baits- this time, use assembly1
#add in a minimum match requirement with -m

#...........................................................................
msg_banner "now assembling nanopore mt reads"

flye --nano-raw $nano_mt2 --genome-size $genome_size \
--out-dir flye002 --threads $threads

assembly2=assembly2.fasta

cp flye002/assembly.fasta $assembly2

#...........................................................................
msg_banner "Round 3"

#...........................................................................
msg_banner "now extracting mitome nanopore reads from all reads"

nano_mt3=nano_mt3.fq.gz

minimap2 -m 5000 -a -x map-ont -t $threads $assembly2 $nano_raw | 
samtools fastq -0 $nano_mt3 -n -F 4 -

#map the reads to a set of baits- this time, use assembly2
#add in a minimum match requirement with -m

#...........................................................................
msg_banner "now assembling nanopore mt reads"

flye --nano-raw $nano_mt3 --genome-size $genome_size \
--out-dir flye003 --threads $threads

assembly3=assembly3.fasta

cp flye003/assembly.fasta $assembly3

#...........................................................................
msg_banner "now polishing flye assembly with long reads"

#round 1. make overlaps file, then run racon

minimap2 -x map-ont -t $threads $assembly3 $nano_mt3 \
| gzip > overlaps1.paf.gz

assembly3_racon = assembly3_racon.fasta

racon --threads $threads $nano_mt3 overlaps1.paf.gz \
$assembly3 > $assembly3_racon

#here, further rounds of polishing made little difference
#option to add in medaka polishing here

#...........................................................................
msg_banner "now extracting illumina mitome reads from all reads"

R1_mt=R1_mt.fq.gz
R2_mt=R2_mt.fq.gz

minimap2 -a -x sr $assembly3_racon $R1 $R2 \
| samtools fastq -1 $R1_mt -2 $R2_mt -F 0x4 -f 0x2 -

#extract cp reads only by mapping to long-read assembly
#needs end dash for stdin
#more info https://broadinstitute.github.io/picard/explain-flags.html
#samtools flags: -F4 exclude unmapped reads, -f2 include properly paired reads

#...........................................................................
msg_banner "now creating subset of illumina mitome reads"

R1_mt_subset = R1_mt_subset.fq.gz
R2_mt_subset = R2_mt_subset.fq.gz

rasusa -i $R1_mt -i $R2_mt --coverage 300 \
--genome-size $genome_size -o $R1_mt_subset -o $R2_mt_subset

#downsample to required coverage, x300

#...........................................................................
msg_banner "now polishing flye-racon assembly with illumina"

#round 1pilon polish

bwa index $assembly3_racon

bwa mem -t $threads $assembly3_racon $R1_mt_subset $R2_mt_subset \
| samtools sort > flye_aln1.bam

samtools index flye_aln1.bam

samtools faidx $assembly3_racon

assembly3_racon_pilon = assembly3_racon_pilon.fasta

pilon --genome $assembly_flye_racon2 --frags flye_aln1.bam \
--output assembly3_racon_pilon \
--fix bases --mindepth 0.5 --changes --threads $threads --verbose

#...........................................................................
msg_banner "now running unicycler assembler"

unicycler -1 $R1_mt_subset -2 $R2_mt_subset -l $nano_mt3 -o unicycler \
--threads $threads --no_rotate --keep 2

#--keep 2 will keep final files but also SAM
#--no_rotate means don't rotate replicons to certain start pos

assembly_unicycler=assembly_unicycler.fasta

cp unicycler/assembly.fasta $assembly_unicycler

#...........................................................................
msg_banner "now getting assembly stats"

seqkit stats $assembly1 $assembly2 $assembly3 \
$assembly3_racon $assembly3_racon_pilon $assembly_unicycler \
-Ta > assembly_stats.tsv

#...........................................................................
msg_banner "Script finished!"
