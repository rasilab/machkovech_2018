# for DNA sequence manipulation
suppressPackageStartupMessages(library(Biostrings))
# for IO fastqfiles
suppressPackageStartupMessages(library(ShortRead))
# for string concatenation
suppressPackageStartupMessages(library(glue))
# for string analysis
suppressPackageStartupMessages(library(stringr))
# for converting between bioconductor and tidyverse
suppressPackageStartupMessages(library(biobroom))
# for munging data
suppressPackageStartupMessages(library(tidyverse))

# get sample name from command line
args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]
print(sample)

# change to scripts directory
setwd("scripts/")

# read in flu genome
flu_genome_file <- "/fh/fast/subramaniam_a/db/rasi/genomes/virus/flu/lowctgnp_wsn_pr8/lowctgnp_wsn_pr8.fasta"
flu_genome <- Biostrings::readDNAStringSet(flu_genome_file) 

# these are the genes to which we match
genes <- c("NPhighCTG", "NPlowCTG", "NP")
# sample <- "ltm_vir"


# read in the trimmed reads for this sample
trimfile <- glue("../processeddata/{sample}/trim.fq")
trimreads <- trimfile %>% 
  FastqFile() %>% 
  readFastq() %>% 
  # extract the sequences as a DNAStringSet
  sread() 
print(trimreads)

# set names for output of read number
names(trimreads) <- seq(length(trimreads))

# get the match locations of aligned reads for each gene 
get_matches <- function(gene){
 matches <- matchPDict(trimreads, flu_genome[[gene]], max.mismatch = 0) %>% 
  unlist()
 GRanges(seqnames = gene, ranges = matches)
}

# write matches for all genes to a tsv file for later analyzing
matches <- map(genes, get_matches) %>% 
  GRangesList() %>% 
  unlist() %>% 
  tidy() %>% 
  select(names, seqname, start, end) %>% 
  write_tsv(glue("../processeddata/{sample}/np_alignments.tsv"))

print("Done alignments")
