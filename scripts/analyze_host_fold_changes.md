Determine list of induced genes under various treatments
================
rasi
02 November, 2018

-   [Import libraries and define analysis specific variables](#import-libraries-and-define-analysis-specific-variables)
-   [Read gene annotations](#read-gene-annotations)
-   [Read in gene counts and process it for input to DESeq2](#read-in-gene-counts-and-process-it-for-input-to-deseq2)
-   [Run DESeq2 to get log fold changes](#run-deseq2-to-get-log-fold-changes)
-   [Get fold-changes between specific sample pairs](#get-fold-changes-between-specific-sample-pairs)
-   [Get list of inferferon-induced genes](#get-list-of-inferferon-induced-genes)
-   [Get list of virus-induced genes](#get-list-of-virus-induced-genes)
-   [Get list of inferferon + virus induced genes](#get-list-of-inferferon-virus-induced-genes)

We look for genes that are induced across both Ribo-seq and RNA-seq samples. Therefore, we treat these two measurements for each treatment as replicates.

Import libraries and define analysis specific variables
=======================================================

``` r
library(DESeq2)
library(tidyverse)

# set the lowest counts to be these for DESeq2 fold-changes if it is zero
deseq2_floor <- 1
```

Read gene annotations
=====================

Select only genes that have annotated protein coding sequences in the CCDS database.

``` r
gene_annotations <- read_tsv(glue::glue("ftp://ftp.ebi.ac.uk/pub/databases/",
  "genenames/new/tsv/locus_groups/protein-coding_gene.txt")) %>%
  filter(!is.na(ccds_id)) %>%
  select(ensembl_gene_id, symbol) %>%
  rename(gene_id = ensembl_gene_id) %>% 
  print()
```

    ## # A tibble: 18,876 x 2
    ##    gene_id         symbol 
    ##    <chr>           <chr>  
    ##  1 ENSG00000121410 A1BG   
    ##  2 ENSG00000148584 A1CF   
    ##  3 ENSG00000175899 A2M    
    ##  4 ENSG00000166535 A2ML1  
    ##  5 ENSG00000184389 A3GALT2
    ##  6 ENSG00000128274 A4GALT 
    ##  7 ENSG00000118017 A4GNT  
    ##  8 ENSG00000094914 AAAS   
    ##  9 ENSG00000081760 AACS   
    ## 10 ENSG00000114771 AADAC  
    ## # ... with 18,866 more rows

Read in gene counts and process it for input to DESeq2
======================================================

Exclude Ribo-seq + LTM data for fold-change calculation

``` r
counts <- list.files("../processeddata/", pattern = "gencode.genes.results",
                          recursive = T,
                          full.names = T) %>% 
  enframe("sno", "file") %>% 
  # exclude LTM data
  filter(!str_detect(file, "ltm")) %>% 
  # get sample name
  mutate(sample = str_extract(file, "[^/]+(?=/gencode.genes.results)")) %>% 
  # read data in
  mutate(data = map(file, data.table::fread)) %>% 
  unnest(data) %>% 
  # remove unwanted columns
  select(-sno, -file, -length, -effective_length, -TPM, -FPKM, -`transcript_id(s)`) %>% 
  mutate(expected_count = as.integer(expected_count)) %>% 
  # first select only gene, sample pairs with non-zero counts
  filter(expected_count > 0) %>% 
  # next create all sample-gene pairs filling with deseq2_floor if it doesn't exist
  complete(sample, gene_id, fill = list("expected_count" = deseq2_floor)) %>% 
  # select only rows that have a minimum of 100 counts across all samples
  group_by(gene_id) %>% 
  mutate(total_counts = sum(expected_count)) %>% 
  ungroup() %>% 
  filter(total_counts > 100) %>% 
  select(-total_counts)

gene_counts <- counts %>% 
  # create one column per sample
  spread(sample, expected_count) %>% 
  # remove gene id suffix to match with gene_name
  mutate(gene_id = str_extract(gene_id, "\\w+"))

print(gene_counts)
```

    ## # A tibble: 11,191 x 9
    ##    gene_id cyclo_ifn cyclo_ifn_vir cyclo_untr cyclo_vir mrna_ifn
    ##    <chr>       <dbl>         <dbl>      <dbl>     <dbl>    <dbl>
    ##  1 ENSG00…       100            90        102        76      315
    ##  2 ENSG00…       263           245        279       175      278
    ##  3 ENSG00…        24             8         14         7       47
    ##  4 ENSG00…        55            54         76        51      151
    ##  5 ENSG00…       323           291         53        43      289
    ##  6 ENSG00…       178           189        175       146      372
    ##  7 ENSG00…       638           567        705       490      903
    ##  8 ENSG00…        41            28         50        20       89
    ##  9 ENSG00…         8            11         11         4       13
    ## 10 ENSG00…         5             2          6         4       27
    ## # ... with 11,181 more rows, and 3 more variables: mrna_ifn_vir <dbl>,
    ## #   mrna_untr <dbl>, mrna_vir <dbl>

``` r
gene_counts <- gene_counts %>% 
  as.data.frame() %>% 
  column_to_rownames("gene_id")
```

Run DESeq2 to get log fold changes
==================================

Treat Ribo and mRNA as replicates to get non-zero P-values

``` r
coldata <- colnames(gene_counts) %>% 
  enframe("sno", "sample") %>% 
  mutate(treatment = str_extract(sample, "(?<=mrna_|cyclo_).+$")) %>% 
  mutate(type = str_extract(sample, "^(mrna|cyclo)")) %>% 
  mutate(colname = sample) %>% 
  as.data.frame() %>% 
  column_to_rownames("colname") %>% 
  select(-sno) %>% 
  print()
```

    ##                      sample treatment  type
    ## cyclo_ifn         cyclo_ifn       ifn cyclo
    ## cyclo_ifn_vir cyclo_ifn_vir   ifn_vir cyclo
    ## cyclo_untr       cyclo_untr      untr cyclo
    ## cyclo_vir         cyclo_vir       vir cyclo
    ## mrna_ifn           mrna_ifn       ifn  mrna
    ## mrna_ifn_vir   mrna_ifn_vir   ifn_vir  mrna
    ## mrna_untr         mrna_untr      untr  mrna
    ## mrna_vir           mrna_vir       vir  mrna

``` r
ddsObject <- DESeqDataSetFromMatrix(countData = gene_counts,
                                    colData = coldata,
                                    design = ~ treatment)

dds <- DESeq(ddsObject)
```

Get fold-changes between specific sample pairs
==============================================

The sample pairs were created manually in `../tables/samplepairs_for_deseq2.tsv`.

``` r
lfc <- read_tsv("../tables/samplepairs_for_deseq2.tsv") %>% 
  mutate(deseq_results = map2(treatment1, treatment2, function(x, y) 
    results(dds, contrast = c("treatment", x, y)))) %>%
  mutate(lfc = map(deseq_results, function(res)
    res %>% as.data.frame() %>% rownames_to_column("gene_id"))) %>%
  select(-deseq_results) %>%
  unnest() %>%
  print()
```

    ## # A tibble: 33,573 x 9
    ##    treatment1 treatment2 gene_id baseMean log2FoldChange lfcSE    stat
    ##    <chr>      <chr>      <chr>      <dbl>          <dbl> <dbl>   <dbl>
    ##  1 ifn        untr       ENSG00…    181.          0.0555 0.570  0.0974
    ##  2 ifn        untr       ENSG00…    267.          0.128  0.741  0.172 
    ##  3 ifn        untr       ENSG00…     24.1         0.318  0.719  0.443 
    ##  4 ifn        untr       ENSG00…     89.9        -0.0123 0.384 -0.0321
    ##  5 ifn        untr       ENSG00…    195.          2.51   0.679  3.69  
    ##  6 ifn        untr       ENSG00…    263.          0.0944 0.290  0.325 
    ##  7 ifn        untr       ENSG00…    756.         -0.0398 0.487 -0.0818
    ##  8 ifn        untr       ENSG00…     52.4        -0.146  0.416 -0.352 
    ##  9 ifn        untr       ENSG00…     11.0        -0.638  0.753 -0.847 
    ## 10 ifn        untr       ENSG00…     12.9        -0.518  1.12  -0.461 
    ## # ... with 33,563 more rows, and 2 more variables: pvalue <dbl>,
    ## #   padj <dbl>

Get list of inferferon-induced genes
====================================

We treat Ribo and mRNA as replicates

``` r
lfc %>% 
  filter(treatment1 == "ifn") %>% 
  left_join(gene_annotations) %>% 
  rename(lfc = log2FoldChange) %>% 
  select(-stat, -lfcSE, -padj) %>% 
  filter(pvalue < 0.001 & lfc >= log2(2)) %>%
  arrange(pvalue) %>% 
  write_tsv("../tables/interferon_induced_genes.tsv") %>% 
  print()
```

    ## # A tibble: 144 x 7
    ##    treatment1 treatment2 gene_id         baseMean   lfc    pvalue symbol
    ##    <chr>      <chr>      <chr>              <dbl> <dbl>     <dbl> <chr> 
    ##  1 ifn        untr       ENSG00000185745    4363.  9.05 3.68e-129 IFIT1 
    ##  2 ifn        untr       ENSG00000119917    1676.  8.21 4.04e- 88 IFIT3 
    ##  3 ifn        untr       ENSG00000135114     620.  7.04 1.22e- 58 OASL  
    ##  4 ifn        untr       ENSG00000156587     411.  4.86 2.47e- 42 UBE2L6
    ##  5 ifn        untr       ENSG00000115267     583.  7.12 2.80e- 36 IFIH1 
    ##  6 ifn        untr       ENSG00000107201    1090.  6.13 3.25e- 28 DDX58 
    ##  7 ifn        untr       ENSG00000168394     543.  3.81 4.60e- 28 TAP1  
    ##  8 ifn        untr       ENSG00000177409     283.  6.10 3.07e- 27 SAMD9L
    ##  9 ifn        untr       ENSG00000119922    1113.  7.37 7.91e- 27 IFIT2 
    ## 10 ifn        untr       ENSG00000138646     510.  4.43 9.27e- 26 HERC5 
    ## # ... with 134 more rows

Get list of virus-induced genes
===============================

We treat Ribo and mRNA as replicates

``` r
lfc %>% 
  filter(treatment1 == "vir") %>% 
  left_join(gene_annotations) %>% 
  rename(lfc = log2FoldChange) %>% 
  select(-stat, -lfcSE, -padj) %>% 
  filter(pvalue < 0.001 & lfc >= log2(2)) %>%
  arrange(pvalue) %>% 
  write_tsv("../tables/virus_induced_genes.tsv") %>% 
  print()
```

    ## # A tibble: 7 x 7
    ##   treatment1 treatment2 gene_id         baseMean   lfc   pvalue symbol
    ##   <chr>      <chr>      <chr>              <dbl> <dbl>    <dbl> <chr> 
    ## 1 vir        untr       ENSG00000185745   4363.   3.49 1.79e-19 IFIT1 
    ## 2 vir        untr       ENSG00000119917   1676.   2.51 1.71e- 8 IFIT3 
    ## 3 vir        untr       ENSG00000213928    149.   2.45 4.18e- 7 IRF9  
    ## 4 vir        untr       ENSG00000135114    620.   2.13 1.01e- 5 OASL  
    ## 5 vir        untr       ENSG00000119922   1113.   2.67 1.67e- 4 IFIT2 
    ## 6 vir        untr       ENSG00000258659     60.8  3.76 2.66e- 4 TRIM34
    ## 7 vir        untr       ENSG00000281394     26.3  7.12 3.59e- 4 <NA>

Get list of inferferon + virus induced genes
============================================

We treat Ribo and mRNA as replicates

``` r
lfc %>% 
  filter(treatment1 == "ifn_vir") %>% 
  left_join(gene_annotations) %>% 
  rename(lfc = log2FoldChange) %>% 
  select(-stat, -lfcSE, -padj) %>% 
  filter(pvalue < 0.001 & lfc >= log2(2)) %>%
  arrange(pvalue) %>% 
  write_tsv("../tables/interferon_plus_virus_induced_genes.tsv") %>% 
  print()
```

    ## # A tibble: 155 x 7
    ##    treatment1 treatment2 gene_id         baseMean   lfc    pvalue symbol
    ##    <chr>      <chr>      <chr>              <dbl> <dbl>     <dbl> <chr> 
    ##  1 ifn_vir    untr       ENSG00000185745    4363.  9.33 3.25e-137 IFIT1 
    ##  2 ifn_vir    untr       ENSG00000119917    1676.  8.64 2.51e- 97 IFIT3 
    ##  3 ifn_vir    untr       ENSG00000135114     620.  7.37 2.72e- 64 OASL  
    ##  4 ifn_vir    untr       ENSG00000156587     411.  4.98 2.08e- 44 UBE2L6
    ##  5 ifn_vir    untr       ENSG00000115267     583.  7.34 2.06e- 38 IFIH1 
    ##  6 ifn_vir    untr       ENSG00000119922    1113.  8.70 8.77e- 37 IFIT2 
    ##  7 ifn_vir    untr       ENSG00000168394     543.  4.02 3.97e- 31 TAP1  
    ##  8 ifn_vir    untr       ENSG00000107201    1090.  6.31 8.80e- 30 DDX58 
    ##  9 ifn_vir    untr       ENSG00000138646     510.  4.77 1.59e- 29 HERC5 
    ## 10 ifn_vir    untr       ENSG00000177409     283.  6.19 5.41e- 28 SAMD9L
    ## # ... with 145 more rows
