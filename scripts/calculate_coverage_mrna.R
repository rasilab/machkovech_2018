# for handling htseq alignments
library(GenomicAlignments)
# for string concatenation
library(glue)
# for string analysis
library(stringr)
# for tab data manipulation
library(tidyverse)
# for converting bioC to tidy objects
library(biobroom)

## 'parallel=TRUE' runs the MAP step in parallel and is currently
## implemented for Unix/Mac only.
register(MulticoreParam(8))

setwd("scripts/")

args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]
reference <- args[2]
coord <- args[3]
inputstrand <- args[4]

# sample <- "mrna_untr"
# reference <- "gencode"
# coord <- "genome"
# inputstrand <- "plus"
outputdir <- "../coverage/"

# set read length dependent trimming from 5' end of read
min_read_length <- 25

# open alignments file
bamfile <- BamFile(glue('../processeddata/{sample}/',
                  '{reference}.{coord}.sorted.bam'))

if (inputstrand == "plus") {
  flag = scanBamFlag(isMinusStrand = F)
} 
if (inputstrand == "minus") {
  flag = scanBamFlag(isMinusStrand = T)
} 
# retrieve alignments along with the RSEM posterior probability
alns <- readGAlignments(bamfile, param = ScanBamParam(tag = 'ZW', flag = flag))
# keep only alignments with ZW > 0
alns <- alns[mcols(alns)$ZW > 0.01 & qwidth(alns) >= min_read_length]

# divide by read length to keep normalized by read counts
mcols(alns)$weight <- mcols(alns)$ZW / qwidth(alns)

# calculate coverage weighted by ZW score from rsem
cvg <- coverage(alns, weight = mcols(alns)$weight)
  
cvg %>% 
  GRanges() %>% 
  tidy() %>% 
  mutate(strand = if_else(inputstrand == "plus", "+", "-")) %>% 
  mutate(score = round(score, 5)) %>% 
  select(-width) %>% 
  write_tsv(glue("{outputdir}/{reference}/",
                 "{sample}.{reference}.{inputstrand}.{coord}.tsv.gz"))

message(glue("exported {sample} {reference} {coord} {inputstrand}.\n"))
