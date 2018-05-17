"""Parses sequences of avian influenza genes.

Sequences are parsed for individual subgroups of birds:
gull, duck, chicken 

Only keeps sequences from strains for which all of these
genes can be found. If there are multiple sequences for a strain
name, takes just one.

These sequences are aligned to the 1918 virus, with
gaps relative to this strain stripped away.

Only keep at most 5 strains per year, randomly chosen.
"""


import re
import random
import pips.fasta
import pips.align


def main():
    """Main body of script."""

    # define variables
    max_per_year = 5 # keep <= this many strains per year
    musclepath = '/Users/jbloom/muscle3.8/'
    infile = 'avian.fasta' # input file of human genes
    outfile = '%s_%s.fasta' # output file
    genes = ['PB2', 'PA', 'NP', 'M1', 'M2', 'NS1', 'NS2']
    bird_types = ['gull', 'duck', 'chicken']
    reffiles = '../1918/1918_%s.fasta'
    gap_cutoff = 0.075 # ignore strains with more than this many gaps
    identity_cutoff = 0.65 # ignore strains with less than this much identity
    excluded_seqs = [ # sequences specified manually for exclusion 
                      'A/chicken/Brescia/1902', # online indicates not from 1902
                      'A/turkey/Minnesota/833/1980', # year discrepancy
                    ]

    # parse for the sequences
    seqs = pips.fasta.Read(infile)
    print "Read a total of %d genes from %s." % (len(seqs), infile)
    headmatch = re.compile('^STRAIN\=(?P<strain>.*?) YEAR\=(?P<year>\d*) SUBTYPE=(?P<subtype>[UunknownMmixed\,HN\-\d]*) SEGMENT\=(?P<segment>[PBA12HNSMF\-]+) ACCESSION\=(?P<accession>.*?) COUNTRY=(?P<country>.*?) HOST\=(?P<host>Avian)$')
    unambiguous_seq_match = re.compile('^[ATGC]+$')
    seqs_d = dict([(bird_type, {}) for bird_type in bird_types])
    i = 0
    for (head, seq) in seqs:
        i += 1
        if not (i % 1000):
            print "Parsing sequence %d of %d." % (i, len(seqs))
        exclude_this_seq = False
        for exclude in excluded_seqs:
            if exclude in head:
                #print "Excluding based on manual specification:\n%s" % head
                exclude_this_seq = True
                continue # manually specified for exclusions
        if exclude_this_seq:
            continue
        seq = seq.upper()
        m = headmatch.search(head)
        if not m:
            print "Failed to match header:\n%s" % head
            continue
        gene = m.group('segment')
        if gene not in genes:
            continue
        if not m.group('year'):
            continue
        year = int(m.group('year'))
        if not unambiguous_seq_match.search(seq):
            continue # has ambiguous nucleotides
        strain = m.group('strain').strip()
        for bird_type in bird_types:
            if bird_type in strain:
                break
        else:
            continue # not one of the types of birds we are considering
        if (strain in seqs_d[bird_type]) and (gene in seqs_d[bird_type][strain]):
            continue # already have this gene for this strain
        refgene = pips.fasta.Read(reffiles % gene)[0]
        a = pips.align.Align([refgene, (head, seq)], musclepath, 'MUSCLE')
        (ident, gaps) = pips.align.PairwiseStatistics(a)
        if ident < identity_cutoff or gaps > gap_cutoff:
            print "Alignment too different, with gaps of %.2f and identities of %.2f for:\n%s" % (gaps, ident, head)
            continue
        a = pips.align.StripGapsToFirstSequence(a)
        try:
            prot = pips.fasta.Translate([a[1]], truncate_incomplete=True, translate_gaps=True) # make sure aligned sequence can be translated
        except ValueError:
            print "Excluding this sequence because it cannot be translated:\n%s" % head
            continue
        if strain in seqs_d[bird_type]:
            seqs_d[bird_type][strain][gene] = (a[1])
            if year != seqs_d[bird_type][strain]['YEAR']:
                raise ValueError("Year mismatch of %d versus %d for strain %s" % (year, seqs_d[bird_type][strain]['YEAR'], strain))
        else: 
            seqs_d[bird_type][strain] = {'YEAR':year, gene:(a[1])}
    for (bird_type, strain_d) in seqs_d.iteritems():
        print "\nFor %s, read genes for %d different strains." % (bird_type, len(strain_d))
    # Now determine the strains with all segments
    for (bird_type, strain_d) in seqs_d.iteritems():
        complete_d = {}
        ncomplete = 0
        for (strain, d) in strain_d.iteritems():
            if len(genes) == len(d) - 1: # complete
                ncomplete += 1
                year = d['YEAR']
                concat_head = 'Concatenated %s STRAIN=%s YEAR=%d' % ('+'.join(genes), strain, year)
                concat_seq = ''.join([d[gene][1] for gene in genes])
                d['concatenated'] = (concat_head, concat_seq)
                if year not in complete_d:
                    complete_d[year] = []
                complete_d[year].append(d)
        print "\nFor %s, found complete gene sets for %d strains from %d years." % (bird_type, ncomplete, len(complete_d))
        if not ncomplete:
            print "Not writing any files for %s because there are no complete gene sets." % (bird_type)
            continue
        sampled_d = {}
        nsampled = 0
        for (year, strainlist) in complete_d.iteritems():
            nsampled += min(len(strainlist), max_per_year)
            if len(strainlist) <= max_per_year:
                sampled_d[year] = strainlist
            else:
                sampled_d[year] = random.sample(strainlist, max_per_year)
        print "\nFor %s, there are %d strains from %d years remaining after sampling just %d per year." % (bird_type, nsampled, len(sampled_d), max_per_year)
        # Write output files
        seqs = dict([(gene, []) for gene in genes + ['concatenated']])
        for (year, strainlist) in sampled_d.iteritems():
            for strain_d in strainlist:
                for gene in genes + ['concatenated']:
                    seqs[gene].append(strain_d[gene])
        for gene in genes + ['concatenated']:
            f = outfile % (bird_type, gene)
            print "Writing %d genes to %s" % (len(seqs[gene]), f)
            pips.fasta.Write(seqs[gene], f)
            n = len(seqs[gene][0][1])
            for seq in seqs[gene]:
                if len(seq[1]) != n:
                    raise ValueError("Length mismatch -- sequences may not be aligned?\n%s\n" % (seq[0], seq[1]))



main() # run the script
