#!/usr/bin/env python
#
#Copyright: This module is copyright 2013 by Lisa Cohen,
#under the BSD License.
#Revision 1.1
#Date: April 15, 2013

"""A module containing bioinformatics functions.

This module is written in Python3 and contains 13 commonly-used bioinformatics functions.
See the program "seq_anal_example.py" for an example. 
These were written in Python 3 and can be used indepently of Biopython.
"""

import re

def seq_length(seq):
    """Prints the length of an oligonucleotide sequence.
    
    """
    seqlength=len(seq)
    print('The length of the oligonucleotide sequence is: ',seqlength)

def seq_count(seq):
    """Counts the number of A, T, G, C bases in a DNA sequence.
    
    The counts are stored and returned in a dictionary with bases as keys and counts as values

    """
    ATGC={'A':0,'T':0,'G':0,'C':0}
    for seq_letter in seq:
           seqletter=seq_letter.upper()
           ATGC[seqletter]+=1     
    return (ATGC)

def GCcont(ATGC):
    """Calculates the GC content of a sequence.
    
    The following equation is used:
    ((G+C)/(A+T+G+C))*100

    """
    GC=((ATGC['G']+ATGC['C'])/(ATGC['A']+ATGC['T']+ATGC['G']+ATGC['C']))*100
    print('The G-C content of this sequence is:',GC,'%')

def ATGCrat(ATGC):
    """Calculates the AT/GC ratio of a sequence. 
    
    The following equation is used:
    (A+T)/(G+C)=AT/GC

    """
    ATGC_ratio=(ATGC['A']+ATGC['T'])/(ATGC['G']+ATGC['C'])
    print('The AT/GC ratio is: ',ATGC_ratio)

def mol_wt(ATGC):
    """Calculates the molecular weight (mw) of a sequence (g/mol).
    
    The calculation uses the mw of each base:
    A: 313.21 g/mol,
    T: 304.2 g/mol,
    C: 289.18 g/mol,
    G: 329.21 g/mol
    The removal of HPO2 and addition of 2H during polymerization is considered in the equation.
    The calculation assumes that the oligonucleotide sequence is more than 8 bases long.
    For more information, please see the complete explanation and oligonucleotide calculator provided by Northwestern University:
    http://www.basic.northwestern.edu/biotools/oligocalc.html

    """
    A_mw=313.21
    T_mw=304.2
    C_mw=289.18
    G_mw=329.21
    mw=(A_mw*ATGC['A'])+(T_mw*ATGC['T'])+(G_mw*ATGC['G'])+(C_mw*ATGC['C'])-61.96
    print('The molecular weight (g/mol) of the sequence is: ',mw)

def Tm(ATGC):
    """Calculates the melting temperature (degC) of the sequence.
    
    The calculation is based on two standard equations taking into account thermodynamic phase transitions of the polymer.
    Standard salt conditions are assumed, with Na+ or K+ between 0.1-1M.
    For more information, please see the complete explanation and oligonucleotide calculator provided by Northwestern University:
    http://www.basic.northwestern.edu/biotools/oligocalc.html

    """
    if sum(ATGC.values())<=14:
        Tempmelt=((ATGC['A']+ATGC['T'])*2)+((ATGC['G']+ATGC['C'])*4)
    else:
        Tempmelt =64.9+(41*(ATGC['G']+ATGC['C']-16.4)/(ATGC['A']+ATGC['T']+ATGC['G']+ATGC['C']))
    print('The melting temperature (degC) under standard salt conditions is: ',Tempmelt)


def transcription(DNA):
    """Transcribes a DNA sequence and returns the mRNA sequence:

    """
    DNA=DNA.upper()
    RNA=DNA.replace('T','U')
    print('The mRNA sequence is:',RNA)
    return RNA
       
def revcomp_RNA(RNA):
    """Reads an RNA sequence and returns its reverse complement:

    """
    revcompRNAseq=RNA.translate(str.maketrans('GCAU','CGUA'))[::-1]
    print('The reverse complement is: ',revcompRNAseq)
    return revcompRNAseq

def revcomp_DNA(DNA):
    """Reads a DNA sequence and returns its reverse complement:

    """
    revcompDNAseq=DNA.translate(str.maketrans('GCAT','CGTA'))[::-1]
    return revcompDNAseq

def codon(RNA):
    """Compiles and returns a dictionary of codons
    
    The RNA sequence is read and stored with the frame (+1, +2, or +3) as keys and codons as values:

    """
    codons={'+1':[],'+2':[],'+3':[]}
    while RNA!='':
        codons['+1'].append(RNA[0:3])
        codons['+2'].append(RNA[1:4])
        codons['+3'].append(RNA[2:5])
        RNA=RNA[3:]
    return codons

def ORF(codons):
    """Open Reading Frame finder (partial):
    
    Finds a partial open reading frame (ORF) or coding sequence (partial cds) 
    by searching for STOP or START (Methionine) codons.
    
    """
    #NOTE: Future versions of this function will piece together a full ORF, if present.
    RNA_ORF={}
    Met='AUG'
    #The codon AUG encodes the amino acid Methionine and acts as an initiation site for translation to START.
    STOP=['UAA','UAG','UGA']
    for strand, codon in codons.items():
        if Met in codon:
            RNA_ORF[strand]=codon[codon.index(Met):]
        for item in STOP:
            if item in codon:
                RNA_ORF[strand]=codon[:codon.index(item)]
    print('These are the strands with open reading frames: ',RNA_ORF)
    return RNA_ORF
    
def aminoacid_codons():
    """Amino Acid Codons: RNA-->protein
    
    
    During the process of translation, ribosomes catalyze the synthesis of polypeptides from amino acids. Individual amino acids are added onto the growing polymer by tRNAs that recognize a specific RNA codon (set of 3 RNA bases).  
    Degeneracy in the genetic code allows for tRNA to recognize several different codons, e.g. GGU, GGC, GGA, and GGG all encode for the amino acid Glycine.
    This function returns a dictionary with amino acid single letter codes (SLC) as keys and RNA codons as values. This dictionary can be used by other functions to look up SLC with RNA codons or vice versa. 
    
    """
    aa={'I':['AUU','AUC','AUA'],'L':['CUU','CUC','CUA','CUG','UUA','UUG'],\
            'V':['GUU','GUC','GUA','GUG'],'F':['UUU','UUC'],'M':['AUG'],\
            'C':['UGU','UGC'],'A':['GCU','GCC','GCA','GCG'],'G':['GGU','GGC','GGA','GGG'],\
            'P':['ACU','ACC','ACA','ACG'],'T':['ACU','ACC','ACA','ACG'],\
            'S':['UCU','UCC','UCA','UCG','ACU','AGC'],'Y':['UAU','UAC'],\
            'W':['UGG'],'Q':['CAA','CAG'],'N':['AAU','AAC'],'H':['CAU','CAC'],\
            'E':['GAA','GAG'],'D':['GAU','GAC'],'K':['AAA','AAG'],\
            'R':['CGU','CGC','CGA','CGG','AGA','AGG'],'X':['UAA','UAG','UGA']}
    #I=Isoleucine
    #L=Leucine
    #V=Valine
    #F=Phenylalanine
    #M=Methionine=START
    #C=Cysteine
    #A=Alanine
    #G=Glycine
    #P=Proline
    #T=Threonine
    #S=Serine
    #Y=Tyrosine
    #W=Tryptophan
    #Q=Glutamine
    #N=Asparagine
    #H=Histidine
    #E=Glutamic acid
    #D=Aspartic acid
    #K=Lysine
    #R=Arginine
    #X=STOP
    return aa

def translation(aa, RNA_ORF):
    """Amino Acid Translation
    
    Function takes in RNA (ORF) and returns a sequence of amino acids. 
    
    """
    
#search aa to find aa[SLC]=codon then append SLClist.append(key) 
    SLCdict={}
    for strand in RNA_ORF:
        SLCdict[strand]="".join([SLC for codon in RNA_ORF[strand] for SLC in aa.keys() if codon in aa[SLC]])
    print('These are the translated amino acid sequences: ',SLCdict)
    return SLCdict


#Future additions to this module will include:
#1. fasta file import
#full ORF prediction
#2. BLAST alignments
#3. Enrichment analysis: input a list of accession, use Fisher's exact test to assign GO terms, output list of significant functional groups represented
#4. Take a fasta file of ESTs, find open reading frame, BLAST, annotate, assign GO terms, output a new fasta file
