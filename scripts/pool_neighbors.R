library(GenomicFeatures)
library(rtracklayer)
library(tidyverse)
library(stringr)
library(RcppRoll)
library(biobroom)

#setwd('scripts/')
args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
# for testing
# file <- "../coverage/gencode/ltm_vir.gencode.plus.genome.bw"

## pools each value in vector with its two adjacent neigbors
pool_neighbors <- function(score_vector) {
  roll_sum(c(0, score_vector, 0), n = 3)
}        

cvg <- read_tsv(file, na = "") %>% 
  # filter with scores greater than 0.01
  filter(score > 0.01) %>% 
  # create a sequence from start to stop for each range
  mutate(pos = map2(start, end, function(x, y) seq(from = x, to = y))) %>% 
  # expand each range to equal its length
  unnest()  %>% 
  # mutate and unnest to create a single pos for each location
  mutate(start = pos, end = pos) %>% 
  select(-pos) %>% 
  GRanges() 

# combine contiguous ranges
reduced_cvg <- GenomicRanges::reduce(cvg, with.revmap = T) 
# assign a group number to each contiguous range
reduced_cvg$group <- seq_along(reduced_cvg)
# map back group number to original coverage along with the length of each group
cvg$group[unlist(reduced_cvg$revmap)] <- rep(reduced_cvg$group, 
                                             lengths(reduced_cvg$revmap))
cvg$group_length[unlist(reduced_cvg$revmap)] <- rep(lengths(reduced_cvg$revmap), 
                                                    lengths(reduced_cvg$revmap))

newfile <- file %>% 
  str_replace("/gencode", "/gencode_pooled") %>% 
  str_replace("/flu", "/flu_pooled") %>% 
  str_replace(".tsv.gz$", ".pooled.tsv.gz")

cvg %>% 
  tidy() %>% 
  group_by(group) %>% 
  mutate(newscore = if_else(group_length > 1, pool_neighbors(score), score)) %>% 
  mutate(newscore = round(newscore, 5)) %>% 
  rename(oldscore = score) %>% 
  rename(score = newscore) %>% 
  write_tsv(newfile) 