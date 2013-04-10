#!/usr/bin/env python

"""
This is an example program for using the module LisaPySeq.
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
LisaPySeq.ORF(LisaPySeq.codon(RNA))
input('Press ENTER to quit')

"""
For test sequences, please see the NCBI Nucleotide database:
http://www.ncbi.nlm.nih.gov/nuccore
Or use the taxonomy broswer for sequences from specific species:
http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Root    

A few examples:
http://www.ncbi.nlm.nih.gov/nuccore/118151389?report=fasta
http://www.ncbi.nlm.nih.gov/nuccore/118343849?report=fasta
http://www.ncbi.nlm.nih.gov/nuccore/18150835?report=fasta
"""

__author__ = "Lisa Cohen"
__maintainer__ = "Lisa Cohen"
__email__ = "lisa.johnson.cohen@gmail.com"
__copyright__ = "Copyright 2013"
__license__ = "BSD"
__version__ = "1.0"
__status__ = "Development"
