#!/usr/bin/python3
# -*- coding: utf-8 -*-
import argparse
import sys


class SeqType(object):
    DNA_IUPAC_letters = "agctbdhkmnrsvwxy"
    DNA_Thymine = "t"
    RNA_Uracile = "u"
    missing_letters = ".nx_- ?\n\b\t\r"

    @classmethod
    def analyse(self, sequence, length_test=50):
        """
        Take a sequence and detect if it contain "dna" or "protein"
        length_test : nb of letters tested to deternimate seq type (default = 50),
        The probability of observing a protein sequence containing only DNA Alphabet in the first twenty residues is almost null
        """
        nb_gap = 0
        rna = False
        typeofseq = "dna"
        sequence = sequence.lower()

        for n, letter in enumerate(sequence):
            if letter in self.missing_letters:
                nb_gap += 1
            else:
                if not (letter in self.DNA_IUPAC_letters):

                    if letter == self.RNA_Uracile:
                        rna = True
                    else:
                        typeofseq = "protein"
                        break

                if rna and (letter == self.DNA_Thymine):
                    rna = False
                    typeofseq = "protein"
                    break

                if n > (length_test + nb_gap):
                    break

        if n + 1 == nb_gap:
            # empty
            return ""

        if rna:
            typeofseq = "rna"

        return typeofseq


def analyse_file(inputfile):
    """
    Take a fasta file and detect if it contain "dna" or "protein"
    """
    typeofseq = ""

    with open(inputfile,newline=None) as input_handle:

        sequence = ""
        first_line = input_handle.readline()
        for line in input_handle:
            if not line.startswith('>') and not line.startswith('#'):
                sequence += line.lower().strip()

            else:
                if sequence or (not line):

                    currentSeqType = SeqType.analyse(sequence)

                    if not currentSeqType:
                        sys.stderr.write("Warning ! Your alignment contains an empty or non-informative sequence\n")
                        currentSeqType = typeofseq

                    elif typeofseq:
                        if typeofseq != currentSeqType:
                            sys.stderr.write("Warning ! Two types of sequences detected\n")
                    typeofseq = currentSeqType
                    sequence = ""

    return typeofseq


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='+', type=str, action="store", default="", help="input fasta file")
    args = parser.parse_args()

    for f in args.file:
        print(analyse_file(f))
