bioinformatics
==============
This repository contains some commonly-used bioinformatics functions written in Python 3. These were mainly written for learning purposes. 

The module [LisaPySeq.py] () can be used without installation of Biopython.

Contact
=======
Email lisa.johnson.cohen@gmail.com or lcohen49@hboi.fau.edu.

Installation instructions
=========================
Download [LisaPySeq.py] () into a working directory. 
   
      import LisaPySeq

An example program is provided in [seq_anal_example.py] ()

Module Contents
===============
The following functions are available in [LisaPySeq.py] ():

1. Prints the length of an oligonucleotide sequence:

         def seq_length(seq):
      
2. Count the number of A, T, G, C bases in a DNA sequence. The counts are stored and returned in a dictionary with bases as keys and counts as values:

         def seq_count(seq):
            ...
            return (ATGC)
         
3. Calculates the GC content of a sequence:

         def GCcont(ATGC):

4. Calculates and prints the AT/GC ratio of a sequence with the following equation: (A+T)/(G+C)

         def ATGCrat(ATGC):
      
5. Calculates and prints the molecular weight (mw) of a sequence (g/mol). The calculation uses the mw of each base (A: 313.21 g/mol, T: 304.2 g/mol, C: 289.18 g/mol, G: 329.21 g/mol) and takes into account the removal of HPO2 and addition of 2H during polymerization. The calculation assumes that the oligonucleotide sequence is more than 8 bases long. For more information, please see the complete explanation and oligonucleotide calculator provided by [Northwestern University][1].  

         def mol_wt(ATGC):
      
6. Calculates and prints the melting temperature (degC) of the sequence under standard salt conditions. The calculation is based on two standard equations taking into account thermodynamic phase transitions of the polymer. Standard salt conditions are assumed, with Na+ or K+ between 0.1-1M. For more information, please see the complete explanation and oligonucleotide calculator provided by [Northwestern University][1].

         def Tm(ATGC):

7. Transcribes a 5'-->3' DNA sequence and returns its mRNA sequence: 

         def transcription(DNA):
            ...
            return RNA

8. Reads an RNA sequence and returns its reverse complement:

         def revcomp_RNA(RNA):
            return revcompRNAseq

9. Reads a DNA sequence and returns its reverse complement:

         def revcomp_DNA(DNA):
            return revcompDNAseq
            

10. Compiles and returns a dictionary of codons from an RNA sequence on with frame (+1, +2, or +3) as keys and codons as values:

         def codon(RNA):
            ...
            return codons

11. Finds a partial open reading frame (ORF) or coding sequence (partial cds) by searching for STOP or START (Methionine) codons. (NOTE: Future versions of this function will piece together a full ORF, if present):

         def ORF(codons):
            ...
            return RNA_ORF

12. Defines amino acid single letter codes (SLC) from RNA codons: RNA-->protein. This function returns a dictionary with SLC as keys and RNA codons as values. This dictionary can be used by other functions to look up SLC with RNA codons or vice versa.

         def aminoacid_codons():
            ...
            return aa
            
13. Translates an RNA sequence (ORF) into a sequence of corresponding amino acids. This function returns a dictionary with SLC sequences for each strand.

         def translation(aa, RNA_ORF):
            ...
            return SLCdict

Credits
=======
These functions are based on several web-based programs we commonly use:
* [Oligonucleotide Properties Calculator] [1]
* [Reverse Complement] [2] 
* [Translation] [3] 
* [Open Reading Frame] [4]
[1]: http://www.basic.northwestern.edu/biotools/oligocalc.html                   "Oligonucleotide Properties Calculator"
[2]: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html     "Reverse Complement"
[3]: http://www.attotron.com/cybertory/analysis/trans.htm                        "Translation"
[4]: http://www.ncbi.nlm.nih.gov/projects/gorf/                                  "Open Reading Frame"

I found this article particularly inspiring:
* [A Quick Guide for Developing Effective Bioinformatics Programming Skills.] [5]

[5]: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000589   "A Quick Guide for Developing Effective Bioinformatics Programming Skills."
