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

sample <- "cyclo_untr"
reference <- "gencode"
coord <- "genome"
inputstrand <- "plus"
outputdir <- "../coverage/"

# set read length dependent trimming from 5' end of read
min_read_length <- 25
max_read_length <- 39
left.trim  <- c(rep(13, 8), rep(14, 7))
names(left.trim) <- as.character(seq(min_read_length, max_read_length))

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
alns <- alns[mcols(alns)$ZW > 0.01]

# length of query
alnLength <- qwidth(alns)
# how much to trim on left side
alnLeftTrim <- left.trim[as.character(alnLength)]
# Width of remaining alignment
alnWidth <- 1
# Subset to alignments that have positive trimmed length
alnGood <- alnWidth > 0 & alnLength >= min_read_length & alnLength <= max_read_length
if (length(alns[alnGood]) > 0) {
  alnStart <- ifelse(
    strand(alns[alnGood]) == "+",
    # trim by left.trim if aln is on + strand
    alnLeftTrim[alnGood] + 1,
    # trim by aln_length - left.trim if aln is on - strand
    alnLength[alnGood] - alnLeftTrim[alnGood])
  psites <- qnarrow(alns[alnGood],
                    start = alnStart,
                    width = alnWidth)
} else {
  # return empty range if no aln pass the above checks
  psites <- GRanges()
}

mcols(psites)$weighted_ZW <- mcols(psites)$ZW / qwidth(psites)

# calculate coverage weighted by ZW score from rsem
cvg <- coverage(psites, weight = mcols(psites)$weighted_ZW)
  
cvg %>% 
  GRanges() %>% 
  tidy() %>% 
  mutate(strand = if_else(inputstrand == "plus", "+", "-")) %>% 
  mutate(score = round(score, 5)) %>% 
  select(-width) %>% 
  write_tsv(glue("{outputdir}/{reference}/",
                 "{sample}.{reference}.{inputstrand}.{coord}.tsv.gz"))

message(glue("exported {sample} {reference} {coord} {inputstrand}.\n"))
