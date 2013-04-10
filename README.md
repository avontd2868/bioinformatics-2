bioinformatics
==============
This repository contains some commonly-used bioinformatics functions in Python. These were mainly written for learning purposes. 
The module "LisaPySeq.py" can be used without installation of Biopython.

Contact
=======
Lisa Cohen
lisa.johnson.cohen@gmail.com
   
Installation instructions
=========================
Download "LisaPySeq.py" into a working directory. Use the "import LisaPySeq" statement. An example program is provided in "seq_anal_example.py"

Module Contents
===============
1. seq_length: returns oligonucleotide length
2. GCcontent
3. ATGCratio
4. mol_wt: returns molecular weight (g/mol)
5. Tm: returns melting temperature under standard salt conditions
6. ORF: returns partial open reading frames
7. revcomp_RNA: returns reverse complement of an RNA sequence
8. revcomp_DNA: returns reverse complement of a DNA sequence
9. transcription: replaces T with U
10. codon: returns codons from an RNA sequence on multiple frames


Credits
=======
These functions are based on several web-based programs commonly used in my work on gene expression profiling:
Oligonucleotide Properties Calculator: http://www.basic.northwestern.edu/biotools/oligocalc.html
Reverse Complement: http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
Translation: http://www.attotron.com/cybertory/analysis/trans.htm
Open Reading Frame: http://www.ncbi.nlm.nih.gov/projects/gorf/
