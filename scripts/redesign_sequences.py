"""Redesigns influenza genes to add / remove alternative start motifs.

First removes CTG, ATG, GTG motifs from all reading frames as best as possible.
These sequences are denoted as "lowCTG" because this is supposed to remove
alternative start sites.

Then adds CTG motifs in reading frame 1 of the PR8-NP-lowCTG sequence to
create PR8-NP-highCTG motif. This sequence should have lots of alternative
starts in the first reading from.

Does all this under the constraint that changes must be synonymous, and that any
new codons that are introduced are chosen to be as frequent as possible in
natural existing sequences, and must exist in at least a threshold number 
of sequences.

Jesse Bloom, 2012. 
Edited by Heather Machkovech, 2018"""


import os
import fasta
import align


def OtherCodons(codon):
    """Given some codon, returns all other codons for the same amino acid."""
    alternatives = []
    assert len(codon) == 3
    aa = fasta.Translate([('codon', codon)])[0][1]
    nts = ['A', 'T', 'C', 'G']
    for nt1 in nts:
        for nt2 in nts:
            for nt3 in nts:
                alternative = "%s%s%s" % (nt1, nt2, nt3)
                if aa == fasta.Translate([('codon', alternative)])[0][1] and codon.upper() != alternative:
                    alternatives.append(alternative)
    return alternatives


def CountCodonOccurrences(seqs, icodon, codons):
    """Counts the number of times a codon occurs in a sequence set.

    seqs -> a list of sequences in (header, sequence) format. These
        sequences should all be of the same length.
    icodon -> the position of a codon found in the sequences of length
        equal to those in seqs, for 1, 2, 3,... numbering of the codons.
    codons -> a list of codons as three nucleotide strings.
    This function counts the number of times that each of the 
        codons in codons occurs at position icodon in the sequences 
        in seqs. The returned variable is a dictionary counts. It is
        keyed by each of the codons in codons, and the values are the
        number of occurrences at position icodon in the sequences in seqs.
    """
    counts = dict([(codon, 0) for codon in codons])
    n = len(seqs[0][1])
    assert n % 3 == 0
    assert 1 <= icodon <= n / 3
    for (head, seq) in seqs:
        assert len(seq) == n
        seqcodon = seq[3 * (icodon - 1) : 3 * icodon]
        for codon in codons:
            if codon.upper() == seqcodon.upper():
                counts[codon] += 1
    return counts


def CountMotifsInFrame(seq, codonmotif, iframe):
    """Counts the number of motifs in a reading frame.

    'seq' is a string giving a nucleotide sequence, with length
        divisible by three.
    'codonmotif' is a 3-nucleotide codon motif.
    'iframe' is a reading frame (1, 2, or 3).
    Returns the number of occurrences of codonmotif in the
        specified reading frame.
    """
    n = 0
    for i in range(len(seq) / 3):
        codon = seq[3 * i + iframe - 1 : 3 * i + 3 + iframe - 1]
        if codon.upper() == codonmotif.upper():
            n += 1
    return n


def RemoveCodonMotif(seq, codonmotif, icodon, iframe, comparison_seqs, count_threshold):
    """Attempts to remove a nucleotide motif of length three (codon length).

    The motif is only removed if it can be done so without changing the protein
        sequence, and by only using other codons that are common in a set
        of comparison sequences. If multiple alternative synonymous codons
        meet the threshold for being sufficiently common, we choose the most
        common one.
    This method assumes that seq[(icodon - 1) * 3 + iframe - 1 : icodon * 3 + iframe - 1]
        == codonmotif
    seq -> a string giving the sequence from which we want to remove the codon.
        It must be a protein-coding sequence, and so have length divisible by 3.
    codonmotif -> a three nucleotide string giving the motif that we want to remove.
    icodon -> the number of the codon in which the motif starts, with
        1 <= icodon <= len(seq) / 3.
    iframe -> the frame in which the motif starts. It can be 1, 2, or 3. For 
        example, if the motif starts at the beginning of codon icodon, then
        iframe would be 1.
    comparison_seqs -> a list of sequences as (header, sequence) 2-tuples.
        Each sequence should be aligned to seq and so be of the same length.
    count_threshold -> the number of times that a codon must be found
        in count_threshold before we would consider making a change to it.
    The returned variable is the original sequence if the motif cannot be removed, or
        a new version of seq in which the motif has been removed.
    """
    seq = seq.upper()
    assert seq[(icodon - 1) * 3 + iframe - 1 : icodon * 3 + iframe - 1] == codonmotif
    c1 = seq[(icodon - 1) * 3 : icodon * 3]
    c1_alts = OtherCodons(c1)
    c1_alts_counts = CountCodonOccurrences(comparison_seqs, icodon, c1_alts)
    c1_alts = [(c1_alts_counts[codon], codon) for codon in c1_alts]
    c1_alts = [(count, codon) for (count, codon) in c1_alts if count >= count_threshold]
    c1_alts.sort()
    c1_alts.reverse() # now holds all acceptable alternatives with most common first
    c1_alts = [tup[1] for tup in c1_alts]
    if iframe == 1:
        # in the same reading frame as the protein
        if c1_alts: # remove motif to most common alternative
            index = (icodon - 1) * 3
            newseq = "%s%s%s" % (seq[ : index], c1_alts[0], seq[index + 3 : ])
        else:
            newseq = seq # motif cannot be removed
    elif iframe in [2, 3]:
        c2 = seq[icodon * 3 : (icodon + 1) * 3]
        c2_alts = OtherCodons(c2)
        c2_alts_counts = CountCodonOccurrences(comparison_seqs, icodon + 1, c2_alts)
        c2_alts = [(c2_alts_counts[codon], codon) for codon in c2_alts]
        c2_alts = [(count, codon) for (count, codon) in c2_alts if count >= count_threshold]
        c2_alts.sort()
        c2_alts.reverse() # now holds all acceptable alternatives with most common first
        c2_alts = [tup[1] for tup in c2_alts]
        index = (icodon - 1) * 3
        for c1_alt in c1_alts:
            alt = "%s%s" % (c1_alt[iframe - 1 : ], c2[ : iframe - 1])
            if alt != codonmotif: # alternative works
                newseq = "%s%s%s%s" % (seq[ : index], c1_alt, c2, seq[index + 6 : ])
                break
        else: # failed to find alternative by changing first codon; trying changing second
            for c2_alt in c2_alts:
                alt = "%s%s" % (c1[iframe - 1 : ], c2_alt[ : iframe - 1])
                if alt != codonmotif: # alternative works
                    newseq = "%s%s%s%s" % (seq[ : index], c1, c2_alt, seq[index + 6 : ])
                    break
            else:
                newseq = seq # motif cannot be removed
    else:
        raise ValueError("Invalid frame of %r." % iframe)
    return newseq # return the new sequence



def main():
    """Main body of the script."""
    # input / output variables
    musclepath = '/Users/jbloom/muscle3.8/'
    outdir = 'redesigned_sequences' # output directory
    remove_targets = { # targets for removing motifs
     # name : (out name, sequence file, comparison file, codonmotifs)
       
        'PR8-NP' : ('PR8-NP-lowCTG', 'sequences/PR8-NP-cDNA.fasta', 'comparison_sequence_sets/NP_aligned.fasta', ['GTG', 'ATG', 'CTG'])
        
    }
    addctg_targets = { # targets for adding CTG motifs in frame 1
     # name : (out name, sequence file, comparison file)
        'PR8-NP-lowCTG' : ('PR8-NP-highCTG', '%s/PR8-NP-lowCTG.fasta' % outdir, 'comparison_sequence_sets/NP_aligned.fasta'),
    }
    count_threshold = 100 # any new codons introduced must be found naturally in at least this many sequences
    spliced_segments = [
      # (spliced name, first gene, its start, its end, second gene, its start, first exon end, second exon start, its end) with all numbering 1, 2, 3 from beginning of spliced gene
      ('PR8-M', 'PR8-M1', 1, 756, 'PR8-M2', 1, 27, 716, 979),
      ('PR8-NS', 'PR8-NS1', 1, 693, 'PR8-NS2', 1, 30, 503, 909),
    ]
    nonoverlap_range = {}  # gives range of non-overlapping reading frame portion
    for tup in spliced_segments:
        nonoverlap_range[tup[1]] = (tup[6], tup[7])
        s2length = tup[8] - tup[7] + tup[6] - tup[5]
        nonoverlap_range[tup[4]] = (max(tup[6], tup[3] - tup[7]), s2length)

    # begin removing motifs
    for (name, (outname, seqfile, comparisonfile, codonmotifs)) in remove_targets.iteritems():
        print "\nRedesigning %s to create %s." % (name, outname)
        originalseq = fasta.Read(seqfile)
        assert len(originalseq) == 1 
        originalseq = originalseq[0]
        (head, seq) = originalseq
        seq = seq.upper()
        assert len(seq) % 3 == 0
        ncodons = len(seq) / 3
        comparison_seqs = fasta.Read(comparisonfile)
        print "Will use %d comparison sequences for the redesign." % len(comparison_seqs)
        for codonmotif in codonmotifs:
            if len(codonmotif) != 3:
                raise IOError("Method currently only works for 3-nucleotide motifs.")
            for icodon in range(1, ncodons + 1):
                for iframe in [3, 2, 1]:
                    if icodon == ncodons and iframe in [2, 3]:
                        continue # don't look past the end of the last codon
                    index = 3 * (icodon - 1) + iframe - 1
                    if seq[index : index + 3] == codonmotif:
                        newseq = RemoveCodonMotif(seq, codonmotif, icodon, iframe, comparison_seqs, count_threshold)
                        assert fasta.Translate([('newseq', newseq)])[0][1] == fasta.Translate([('seq', seq)])[0][1], "icodon = %d, iframe = %d, seq[index : index + 3] = %s, newseq[index : index + 3] = %s, seq[3 * (icodon - 1) : 3 * (icodon + 1)] = %s, newseq[3 * (icodon - 1) : 3 * (icodon + 1)] = %s" % (icodon, iframe, seq[index : index + 3], newseq[index : index + 3], seq[3 * (icodon - 1) : 3 * (icodon + 1)], newseq[3 * (icodon - 1) : 3 * (icodon + 1)]) 
                        seq = newseq
        # if the gene has overlapping reading frames, only use recoded portions of non-overlapping
        if name in nonoverlap_range:
            print "Adjusting redesign of %s to account for the overlapping reading frames." % name
            (start_no, end_no) = nonoverlap_range[name]
            (startcodon_no, endcodon_no) = (start_no / 3 + 1, end_no / 3 - 1)
            seq = "%s%s%s" % (originalseq[1][ : 3 * startcodon_no], seq[3 * startcodon_no : 3 * endcodon_no], originalseq[1][3 * endcodon_no : ])
        assert fasta.Translate([originalseq])[0][1] == fasta.Translate([('seq', seq)])[0][1] # make sure the protein sequence is unchanged
        header = "%s: %s redesigned to eliminate possible alternative start motifs. The protein sequence is unchanged, and codons are only used if they occur in %d of %d natural sequences." % (outname, name, count_threshold, len(comparison_seqs))
        for codonmotif in codonmotifs:
            for iframe in [1, 2, 3]:
                s = "%s in frame %d reduced from %d to %d." % (codonmotif, iframe, CountMotifsInFrame(originalseq[1], codonmotif, iframe), CountMotifsInFrame(seq, codonmotif, iframe))
                print s
                header = "%s %s" % (header, s)
        outfile = '%s/%s.fasta' % (outdir, outname)
        print "Writing %s to %s." % (outname, outfile)
        fasta.Write([(header, seq)], outfile)

    # now add CTG motifs
    targetcodon = 'CTG'
    targetaa = fasta.Translate([('codon', targetcodon)])[0][1]
    for (name, (outname, seqfile, comparisonfile)) in addctg_targets.iteritems():
        print "\nRedesigning %s to create %s." % (name, outname)
        originalseq = fasta.Read(seqfile)
        assert len(originalseq) == 1 
        originalseq = originalseq[0]
        (head, seq) = originalseq
        seq = seq.upper()
        assert len(seq) % 3 == 0
        ncodons = len(seq) / 3
        comparison_seqs = fasta.Read(comparisonfile)
        print "Will use %d comparison sequences for the redesign." % len(comparison_seqs)
        for icodon in range(ncodons):
            codon = seq[icodon * 3 : icodon * 3 + 3]
            aa = fasta.Translate([('codon', codon)])[0][1]
            if aa == targetaa and codon != targetcodon: # a replacement possibility
                ctgcounts = CountCodonOccurrences(comparison_seqs, icodon + 1, [targetcodon])[targetcodon]
                if ctgcounts >= count_threshold: # replacement possible
                    newseq = "%s%s%s" % (seq[ : icodon * 3], targetcodon, seq[icodon * 3 + 3 : ])
                    assert fasta.Translate([('newseq', newseq)])[0][1] == fasta.Translate([('seq', seq)])[0][1] # make sure the protein sequence is unchanged
                    seq = newseq
        assert fasta.Translate([originalseq])[0][1] == fasta.Translate([('seq', seq)])[0][1] # make sure the protein sequence is unchanged
        s = "%s in frame 1 increased from %d to %d." % (targetcodon, CountMotifsInFrame(originalseq[1], codonmotif, 1), CountMotifsInFrame(seq, codonmotif, 1))
        print s
        header = "%s: %s redesigned to add %s motifs in frame 1. The protein sequence is unchanged, and codons are only used if they occur in %d of %d natural sequences." % (outname, name, targetcodon, count_threshold, len(comparison_seqs))
        header = "%s %s" % (header, s)
        outfile = '%s/%s.fasta' % (outdir, outname)
        print "Writing %s to %s." % (outname, outfile)
        fasta.Write([(header, seq)], outfile)
    # write full alternatively spliced segment coding regions
    for (name, g1, g1_start, g1_end, g2, g2e1_start, g2e1_end, g2e2_start, g2e2_end) in spliced_segments:
        outname = "%s-lowCTG" % name
        (g1name, g2name) = ('%s-lowCTG' % g1, '%s-lowCTG' % g2)        
        g1_file = '%s/%s.fasta' % (outdir, g1name)
        g2_file = '%s/%s.fasta' % (outdir, g2name)
        outfile = "%s/%s.fasta" % (outdir, outname)
        print "Creating %s from %s and %s." % (outname, g1name, g2name)
        (g1_head, g1_seq) = fasta.Read(g1_file)[0]
        (g2_head, g2_seq) = fasta.Read(g2_file)[0]
        merged = "%s%s" % (g1_seq[ : g2e2_start - 1], g2_seq[g2e1_end : ])
        assert g1_seq in merged
        a = align.Align([(outname, merged), (g1name, g1_seq), (g2name, g2_seq)], musclepath, 'MUSCLE')
        a = align.AddDots(a)
        print "Here is the alignment:\n>%s\n%s\n>%s\n%s\n>%s\n%s" % (a[0][0], a[0][1], a[1][0], a[1][1], a[2][0], a[2][1])
        print "Writing %s to %s" % (outname, outfile)
        header = "%s, made by merging the following: %s; %s" % (outname, g1_head, g2_head)
        fasta.Write([(header, merged)], outfile)



main() # run the script
