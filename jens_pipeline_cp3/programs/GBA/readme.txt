Instruction for running GBA
----------------------------

1. Download lcr.zip
2. Unzip it.
3. cd lcr
4. Use the command  java applications.Gbm <protein sequence filename>  <forget rate file> <t2> <t1> <t3>  <substitution matrix file> <output format> 2>/dev/null

Note:
1. <executable>: The main executable is Gbm.class under the package 'applications'. So call applications.Gbm from lcr folder.

2. <protein sequence filename>: Write the protein sequences in fasta format in a file. You can provide more than one proteins in that file.
There needs to be a blank line between two proteins.

3. The <forget rate file> contains the relevant tabular data for a particular data. We keep the useful forget rate file in the folder
swissprotLearnedMatrices/wForgetRate/normalized/

Use the following pathname s for <forget rate file> .

a) For forget rate 0.95
     swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095
b) For forget rate 0.90
    swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095
c) For forget rate 0.85
    swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow085

4. <t1>, <t2> and <t3> are all positive integer values. For more information please see
        http://bioinformatics.cise.ufl.edu/GBA/lcr/Help.htm

5. <substitution matrix file>: This contains an amino acid substitution matrix file such as
Blosum 62. For our purpose we use the file knowledge/blosum62Matrix relative to 'lcr' directory

6. <output option>: For 0 it will provide the range of the repeats in the protein sequence. For 1 it will provide protein
sequence with the repeats masked with 'x'.

7. To avoid the warning message use 2>/dev/null at the end.

So, a sample query from the lcr directory is

java applications.Gbm example/sampleSeq  swissprotLearnedMatrices/wForgetRate/normalized/combinedMatricesRowByRow095 3 15 5  knowledge/blosum62Matrix 1 2>/dev/null

Output
------------
>PER_DROTE
EGSGGSGSSGNFTTASNIHMSSVTNTSIAxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx 
