# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Benjamin Ziemann Github: zneb97

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

    Rationale: ensuring basic return of nucleotides for A/T
    >>> get_complement('A')
    'T'

    Rationale: ensuring basic return of nucleotides
    >>> get_complement('C')
    'G'


    """
    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'T':
        return 'A'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'

    pass


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna = dna[::-1]
    parse = list(dna)
    reverse = ""
    for i in range(len(dna)):
        reverse = reverse + get_complement(parse[i])
    return reverse

    pass


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """

    #Break dna string into triplets
    triples = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    breakpoint = 0

    #Search for stop codons
    for i in range(len(triples)):
        if(triples[i] == "TAG") or (triples[i] == "TGA") or (triples[i] == "TAA"):
            breakpoint = i
            break

    #Return appropriate whole string
    if(breakpoint == 0):
        return dna
    else:
        dna = dna[0:breakpoint*3]
        return dna

    pass


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    starts = []
    orfs = []
    count = 0
    #Break into triples
    triples = [dna[i:i + 3] for i in range(0, len(dna), 3)]
    for i in range(len(triples)):
        if triples[i] == "ATG":
            #If first start found
            if count == 0:
                starts.append(i)
                count = 1
            #Ensure no nesting
            else:
                for q in range(0,len(triples[0:i])):
                    if(triples[q] == "TAG") or (triples[q] == "TGA") or (triples[q] == "TAA"):
                        starts.append(i)
                        break

    #Create list of full orfs
    for j in range(len(starts)):
        orf = dna[starts[j]*3:len(dna)]
        orfs.append(rest_of_ORF(orf))
    return orfs
    pass


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    frameOne = find_all_ORFs_oneframe(dna)
    dna = dna[1:len(dna)]
    frameTwo = find_all_ORFs_oneframe(dna)
    dna = dna[1:len(dna)]
    frameThree = find_all_ORFs_oneframe(dna)

    return frameOne + frameTwo + frameThree
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverseDNA = get_reverse_complement(dna)
    reverseORFs = find_all_ORFs(reverseDNA)
    dnaORFs = find_all_ORFs(dna)
    orfs = dnaORFs + reverseORFs

    return orfs

    pass


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

if __name__ == "__main__":
    import doctest
    #doctest.run_docstring_examples(find_all_ORFs_both_strands, globals(), verbose=True)
    doctest.testmod(verbose=True)
