BEGFE: Bayesian Estimation for Gene Family Evolution

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program implements a Markov Chain Monte Carlo algorithm to estimate the posterior probability distribution of the birth and death rate parameter and the numbers of gene copies at the internodes of the phylogenetic tree.  In addition, BEGFE can simulate gene family data under the birth and death model.

Compiling the program 
To compile the program from source code, type make and hit return under the directory src. 

Running the program
./begfe controlfile

Examples: 
Two control files are included in the package. The control file controlsim is used to simulate gene family data. 
./begfe controlsim

This will produce two files; sim1 and sim1.true. The simulated gene family dataset is saved in the file sim1, while sim1.true contains the true values of the parameters in the Bayesian model. 
The other control file control1 is for carrying out the Bayesian analysis of the gene family data. Type ./begfe control1 and hit return. 

Output
There are two output files; sim1.out and sim1.pvalue. The MCMC output for all parameters is saved in sim1.out. Each column in sim1.out represents the posterior distribution of a particular parameter in the birth and death model. The parameters such as the birth and death rate are estimated by the Bayesian means, i.e., the averages of the columns after discarding the burn-in period. 

The Bayesian p-values (PPP) for all gene families in the dataset are saved in sim1.pvalue. The average of a column in sim1.pvalue is the Bayesian p-values for a particular gene family. Of course, the burn-in period must be discarded.


Citation: Liu, L., L. Yu, V. Kalavacharla, Z. Liu. A Bayesian model for gene family evolution. BMC Bioinformatics. 2011, 12:426
