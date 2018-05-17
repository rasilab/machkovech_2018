"""Parses sequences of human influenza genes.

Gets the PB2, PA, NP, M1, M2, NS1, and NS2 coding
sequences descended from the 1918 virus. Builds these histories
using H1N1 from 1918 to 1957, H2N2 from 1957 to 1968, and
H3N2 from 1968 forward.

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
    infile = 'human.fasta' # input file of human genes
    outfile = 'human_%s.fasta' # output file
    genes = ['PB2', 'PA', 'NP', 'M1', 'M2', 'NS1', 'NS2']
    subtypes = { # keyed by subtype, values are (startyear, endyear)
            'H1N1':(1918, 1957),
            'H2N2':(1957, 1968),
            'H3N2':(1968, 2012),
               }
    reffiles = '../1918/1918_%s.fasta'
    gap_cutoff = 0.075 # ignore strains with more than this many gaps
    identity_cutoff = 0.85 # ignore strains with less than this much identity
    excluded_seqs = [ # sequences specified manually for exclusion after building
                      # phylogenetic trees -- they appear to be mislabeled
                      # or reassortants
                      'A/Indiana/08/2011', 
                      'A/Iowa/08/2011', 
                      'A/Indiana/10/2011', 
                      'A/Melbourne/1/1946', 
                      'A/Nanjing/49/1977', 
                      'A/Moscow/1019/1965', 
                      'A/Ontario/RV1273/2005', 
                      'A/Maine/06/2011', 
                      'A/West Virginia/06/2011', 
                      'A/Indiana/08/2011', 
                      'A/Iowa/07/2011', 
                      'A/Iowa/09/2011', 
                      'A/Pennsylvania/11/2011', 
                      'A/Pennsylvania/09/2011', 
                      'A/Leningrad/1954/1', 
                      'A/Hong Kong/1774/99', 
                      'A/Fujian/411/02-like',
                      'A/WSN/1933 TS61', # redundant with A/WSN/1933
                    ]

    # parse for the sequences
    seqs = pips.fasta.Read(infile)
    print "Read a total of %d genes from %s." % (len(seqs), infile)
    headmatch = re.compile('^STRAIN\=(?P<strain>.*?) YEAR\=(?P<year>\d*) SUBTYPE=(?P<subtype>[UunknownMmixed\,HN\-\d]+) SEGMENT\=(?P<segment>[PBA12HNSMF\-]+) ACCESSION\=(?P<accession>.*?) COUNTRY=(?P<country>.*?) HOST\=(?P<host>Human)$')
    unambiguous_seq_match = re.compile('^[ATGC]+$')
    strain_d = {}
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
            raise ValueError("Failed to match header:\n%s" % head)
        gene = m.group('segment')
        if gene not in genes:
            continue
        subtype = m.group('subtype')
        if subtype not in subtypes:
            continue
        year = int(m.group('year'))
        if not (subtypes[subtype][0] <= year <= subtypes[subtype][1]):
            continue
        if not unambiguous_seq_match.search(seq):
            continue # has ambiguous nucleotides
        strain = m.group('strain').strip()
        if (strain in strain_d) and (gene in strain_d[strain]):
            continue
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
        if strain in strain_d:
            strain_d[strain][gene] = (a[1])
            if year != strain_d[strain]['YEAR']:
                raise ValueError("Year mismatch of %d versus %d for strain %s" % (year, strain_d[strain]['YEAR'], strain))
        else: 
            strain_d[strain] = {'YEAR':year, gene:(a[1])}
    print "\nRead genes for a total of %d different strains." % len(strain_d)
    # Now determine the strains with all segments
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
    print "\nFound complete set of genes for %d different strains from %d years." % (ncomplete, len(complete_d))
    sampled_d = {}
    nsampled = 0
    for (year, strainlist) in complete_d.iteritems():
        nsampled += min(len(strainlist), max_per_year)
        if len(strainlist) <= max_per_year:
            sampled_d[year] = strainlist
        else:
            sampled_d[year] = random.sample(strainlist, max_per_year)
    print "\nThere are %d strains from %d years remaining after sampling just %d per year." % (nsampled, len(sampled_d), max_per_year)
    # Write output files
    seqs = dict([(gene, []) for gene in genes + ['concatenated']])
    for (year, strainlist) in sampled_d.iteritems():
        for strain_d in strainlist:
            for gene in genes + ['concatenated']:
                seqs[gene].append(strain_d[gene])
    for gene in genes + ['concatenated']:
        f = outfile % gene
        print "Writing %d genes to %s" % (len(seqs[gene]), f)
        pips.fasta.Write(seqs[gene], f)
        n = len(seqs[gene][0][1])
        for seq in seqs[gene]:
            if len(seq[1]) != n:
                raise ValueError("Length mismatch -- sequences may not be aligned?\n%s\n" % (seq[0], seq[1]))



main() # run the script
