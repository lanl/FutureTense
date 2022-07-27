# FutureTense
Probabilistic prediction of viral evolution

Version 0.98

## Introduction
The FutureTense software is a collection of scripts and compiled programs for predicting viral genome sequence evolution. Before any predictions can be made (and tested), genome sequence data must be downloaded and preprocessed to infer a genealogy tree representing the hierarchy of parent and child genomes. A reference genome is also required to define the codon positions of each nucleotide in the genome. Once a genealogy tree has been computed, the `motif` program is used to train (and validate) a probabilistic model of genome sequence evolution. Specifically, each allowed mutation (currently single nucleotide substitution, single nucleotide deletion and insertion of one or more nucleotides) at each position in a parent genome is assigned a score that is used to rank the likelihood that a potential mutation will be observed in a child genome sequence. The accuracy of this ranking is assessed by computing the Mann-Whitney U statistic (i.e., Area Under the Reciever Operating Characteritic curve) of child genomes that have been held out of the training set.

### Methodology

## Data download and preprocessing

## Compiling the program
