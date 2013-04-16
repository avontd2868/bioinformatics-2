#!/usr/bin/env python
#
#Copyright: This module is copyright 2013 by Lisa Cohen,
#under the BSD License
#Revision 1.1
#Date: April 15, 2013

"""
Example program using LisaPySeq

This program imports and refers to functions in the LisaPySeq.py module. A user-defined nucleotide sequence is required.
"""

import re, LisaPySeq


seq=input('Please enter any DNA sequence: ')
LisaPySeq.seq_length(seq)
seqcount=LisaPySeq.seq_count(seq)
LisaPySeq.GCcont(seqcount)
LisaPySeq.ATGCrat(seqcount)
LisaPySeq.mol_wt(seqcount)
LisaPySeq.Tm(seqcount)
RNA=LisaPySeq.transcription(seq)
LisaPySeq.revcomp_RNA(RNA)
RNA_ORF=LisaPySeq.ORF(LisaPySeq.codon(RNA))
LisaPySeq.translation(LisaPySeq.aminoacid_codons(),RNA_ORF)
input('Press ENTER to quit')

#For test sequences, please see the NCBI Nucleotide database:
#http://www.ncbi.nlm.nih.gov/nuccore
#Or use the taxonomy broswer for sequences from specific species:
#http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Root    
#
#A few examples:
#http://www.ncbi.nlm.nih.gov/nuccore/118151389?report=fasta
#http://www.ncbi.nlm.nih.gov/nuccore/118343849?report=fasta
#http://www.ncbi.nlm.nih.gov/nuccore/18150835?report=fasta
#
#NOTE: Future versions of this program will input fasta files
#for functions from the LisaPySeq.py to process
