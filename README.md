# BEGFE2.0: Bayesian Estimation for Gene Family Evolution

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program implements a Markov Chain Monte Carlo algorithm to estimate the posterior probability distribution of the birth and death rate parameter and the numbers of gene copies at the internodes of the phylogenetic tree.  In addition, BEGFE can simulate gene family data under the birth and death model.

## New features in version 2.0
1. Parallel computing using PTHREAD is available in version 2.0. Set PTHREAD=yes in Makefile
2. The constraints on the birth rate parameter lambda in the old versions are removed in version 2.0
3. The birth-death probability is calculated using two formulas - one for alpha < 0.5 and another for alpha > 0.5 where alpha = lambda\*brlens / (1+lambda\*brlens)
4. The program provides a warning message if the species names in the data file do not match the names in the control file
5. The program provides a warning message if the number of gene families in the control file > the number of gene families in the data file
6. The second column "NUMBER" is removed in the data file


## Compile from source code 
To compile the program from source code, type make and hit return under the directory src. 

## Run the program
./begfe controlfile

## Example control files 
Two control files are included in the package. The control file controlsim is used to simulate gene family data. 

### A control file for simulating gene family data
1 #0:analysis, 1:simulation
sim1 #output file
-1 #seed, -1:random seed
124 #number of gene families
5 #number of species
(((chimp:6#0.005,human:6#0.005):81#0.005,(mouse:17#0.005,rat:17#0.005):70#0.005):6#0.005,dog:93#0.005)#0.005; #species tree, the numbers after '#' are birth/death rates
5 15 #the number of gene copies at the tree root is simulated from an uniform distribution (5, 15)

To run the simulation, type ./begfe controlsim

This will produce two files; sim1 and sim1.true. The simulated gene family dataset is saved in the file sim1, while sim1.true contains the true values of the parameters in the Bayesian model. 

### A control file for analyzing gene family data

The other control file control1 is for carrying out the Bayesian analysis of the gene family data. Type ./begfe control and hit return. 

0 #0:analysis, 1:simulation
sim1 #input file
-1 #seed, -1:random seed
124 #number of gene families
5 #number of species
(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:93); #species tree
10000 100 0 #number of generations, save every 100 trees, 0:unlink and 1:link

There are two output files; sim1.out and sim1.pvalue. The MCMC output for all parameters is saved in sim1.out. Each column in sim1.out represents the posterior distribution of a particular parameter in the birth and death model. The parameters such as the birth and death rate are estimated by the Bayesian means, i.e., the averages of the columns after discarding the burn-in period. 

The Bayesian p-values (PPP) for all gene families in the dataset are saved in sim1.pvalue. The average of a column in sim1.pvalue is the Bayesian p-values for a particular gene family. Of course, the burn-in period must be discarded.


## Citation 
Liu, L., L. Yu, V. Kalavacharla, Z. Liu. A Bayesian model for gene family evolution. BMC Bioinformatics. 2011, 12:426

## Old versions
Old versions BEGFE are available at https://github.com/lliu1871/oldversion

## Example data set
The Anolis MHC example data set (unpublished as of April 2022) is provided courtesy of Daren Card and Scott Edwards, Harvard University


