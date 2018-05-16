#' Combine fastq files and rename them to a nice format

library(ShortRead)
library(stringr)
library(glue)
library(tidyverse)

args <- commandArgs(trailingOnly = T)
sampleindex <- as.integer(args[1])
# sampleindex <- 1

rawfiles <- list.files("../rawdata/", 
                       pattern = ".+fastq.gz$",
                       full.names = T, recursive = T)
samples <- str_extract(rawfiles, "[^/]+(?=_[ACTG]{6})") %>% 
  unique()

sample <- samples[sampleindex]
outputfile <- str_replace_all(sample, "-", "_")
outputfile <- glue("../rawdata/{outputfile}.fq.gz")

rawfiles <- rawfiles %>% str_subset(paste0(sample,"_[ACTG]{6}")) %>% FastqStreamerList(n = 1e6)


writetofile <- function(file) {
  stream <- open(file)
  on.exit(close(stream))
  
  repeat{
    fq <- yield(stream)
    if (length(fq) == 0)
      break
    writeFastq(fq, outputfile, "a", compress = T)
  }
}

temp <- lapply(rawfiles, writetofile)
print(rawfiles)
