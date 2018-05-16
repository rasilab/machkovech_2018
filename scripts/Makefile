# -*- Mode: makefile -*-
SHELL=/bin/bash -x

SAMPLE=ENTER_SAMPLE_FROM_COMMAND_LINE
OUTPUTDIR:=processeddata/$(SAMPLE)

trim:  $(OUTPUTDIR)/trim.fq
removerrna: $(OUTPUTDIR)/trim.fq
rsem_gencode: $(OUTPUTDIR)/gencode.genome.sorted.bam
rsem_flu: $(OUTPUTDIR)/flu.genome.sorted.bam

PYTHON=/fh/fast/subramaniam_a/user/rasi/virtualenv/default2/bin/python

FASTX_RC=/fh/fast/subramaniam_a/user/rasi/bin/biotools/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/fastx_reverse_complement

CUTADAPT=/fh/fast/subramaniam_a/user/rasi/virtualenv/default3/bin/cutadapt

BOWTIE_PATH="/app/bowtie/1.1.1"
BOWTIE=$(BOWTIE_PATH)/bowtie

RSEM_ALIGN=/fh/fast/subramaniam_a/user/rasi/bin/rsem/rsem-calculate-expression

RSEM_WIGGLE=/fh/fast/subramaniam_a/user/rasi/bin/rsem/rsem-bam2wig

RSCRIPT=/fh/fast/subramaniam_a/user/rasi/src/R-3.3.2/bin/Rscript

RRNA_BOWTIE_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/hg38.rrna
RSEM_GENCODE_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/hg38.gencode.v24.rsem/bowtie
RSEM_FLU_INDEX=/fh/fast/subramaniam_a/db/rasi/bowtie/flu.pr8.wsn.recodednp.rsem/bowtie

# reads are trimmed to a final length of 40 after a poly A adapter of 10 nt or longer is removed
CUTADAPT_ARGS=--adapter=AAAAAAAAAA --length=40 --minimum-length=25 

# for ribosome profiling samples, reads without adapters or outside the 25-35nt window are discarded 
CUTADAPT_RIBO_EXTRA_ARGS=--discard-untrimmed --maximum-length=40 

# the q33 is required for standard illumina fastq output files
FASTX_RC_ARGS=-Q33

RRNA_BOWTIE_ARGS=--seedlen=23 --threads=8

RSEM_ARGS=--num-threads 8 --output-genome-bam --sort-bam-by-coordinate  --quiet-em
# make seed length shorter to prevent short reads from being discarded, allow  3 mismatches in seed region
RSEM_FLU_EXTRA_ARGS=--seed-length 21 --bowtie-n 3

## all: run analysis to generate all final files
all: trim
all: removerrna 
all: rsem_gencode
all: rsem_flu

## removetrim: removes trimmed fastq files
removetrim: 
	rm --force $(OUTPUTDIR)/trim.fq

## removerrna: remove fq files depleted of rrna contaminants
removenorrna:
	rm --force $(OUTPUTDIR)/norrna.fq

## removersem: remove rsem output files
removersem:
	rm -r --force processeddata/*/gencode.transcripts*bam

## clean: removes all processed data 
clean: removetrim
clean: removerrna
clean: removersem

## help: displays this help
help: Makefile
	@sed -n 's/^##//p' $<

# prevents error in case files are named with these keywords
.PHONY: clean all help

# trim polyA adapters (for RP reads) 
# additionally, reverse complement for mrna reads before trimming
processeddata/%/trim.fq: rawdata/%.fq.gz
	mkdir -p processeddata/$*
# for mrna samples
ifneq (,$(findstring mrna, $(SAMPLE)))
	gunzip -c $< | \
	  	$(FASTX_RC) $(FASTX_RC_ARGS) | \
		$(CUTADAPT) $(CUTADAPT_ARGS) --output $@ - \
		1> $@.log 2> $@.log
# for ribosome profiling samples
else
	$(CUTADAPT) $(CUTADAPT_ARGS) $(CUTADAPT_RIBO_EXTRA_ARGS) --output $@ $< \
		1> $@.log 2> $@.log
endif

# align against rrna contaminants
processeddata/%/norrna.fq: processeddata/%/trim.fq
	$(BOWTIE) \
	$(RRNA_BOWTIE_ARGS) --un $@ --sam $(RRNA_BOWTIE_INDEX) $< \
	2> $@.log \
	> /dev/null

# align norrna reads against gencode transcripts using rsem
%/gencode.genome.sorted.bam: %/norrna.fq
	$(RSEM_ALIGN) \
	$(RSEM_ARGS) --bowtie-path $(BOWTIE_PATH) $< $(RSEM_GENCODE_INDEX) $*/gencode \
	1> $*/gencode.log \
	2> $*/gencode.log
	
# align norrna reads against flu transcripts using rsem
%/flu.genome.sorted.bam: %/norrna.fq
	$(RSEM_ALIGN) \
	$(RSEM_ARGS) $(RSEM_FLU_EXTRA_ARGS) --bowtie-path $(BOWTIE_PATH) $< $(RSEM_FLU_INDEX) $*/flu \
	1> $*/flu.log \
	2> $*/flu.log
