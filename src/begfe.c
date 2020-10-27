/*
 *  genedu 1.1
 *
 *  copyright 2009-2011
 *
 *  Liang Liu
 *  Department of Organismic and Evolutionary Biology
 *  Harvard University
 *
 *  lliu@oeb.harvard.edu
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details (www.gnu.org).
 *
 * 
 */

#include	"begfe.h"
#include	"tool.h"
#               if defined (MPI_ENABLED)
#include        <mpi.h>
#               endif

void 		InitialParam (Tree *tree);
double		LikelihoodBD (int s, int c, double lambda, double brlens);
int		    Loglike1tree (Tree *tree, int igene, double *loglike);
int 		LoglikeSimtree (Tree *tree, double *loglike);
int 		McmcRun (Tree *tree);
int 		MoveLambda (Tree *tree);
int 		MoveNode (Tree *tree);
void 		PrintHeader (void);
int		    PrintState (int round, FILE *outfile);
int		    ReadData (FILE *fin);
int         Reject (int *numreject);
int 		SimGeneFamily (int inode, int ncopies, Tree *tree, double lambda);
int 		SimNcopies (int s, double lambda, double brlens);
void 		SimNcopiesProb (int s, double lambda, double brlens, int n, double *logprob);
int			Simulation (Tree *tree);

/*local variables*/
Tree 		sptree;
long int	seed=0;
long int	ngenefamily;
double		curLn;
FILE		*fout;
FILE        *fpvalue;
FILE		*fsim;
Chain	 	mcmc;
int			lambdalink;
int			sim;
int			simngene[2];
double		*simlambda;

#       if defined (MPI_ENABLED)
int             proc_id;
int             num_procs;
#       endif

int main (int argc, char *argv[])
{
	FILE *fin;
	clock_t		previousCPUTime, currentCPUTime;

	
	previousCPUTime = clock();

#	if defined (MPI_ENABLED)
		MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);
        MPI_Comm_size(MPI_COMM_WORLD,&num_procs);
	if( proc_id == 0) 
		PrintHeader();
#	else
        PrintHeader();
#	endif

	fin = (FILE*)gfopen(argv[1],"r");
	if(ReadData (fin) == ERROR)
	{
		MrBayesPrint ("Problem in ReadData\n");
		exit (-1);
	}
	fclose (fin);	

	if(sim == 1) {
		if(Simulation (&sptree) == ERROR){
			printf("Errors in MCMCRUN\n");
			exit(-1);
		}
	}else{
		if(McmcRun (&sptree) == ERROR){
			printf("Errors in MCMCRUN\n");
			exit(-1);
		}
	}
	
	/*free memory*/

#	if defined (MPI_ENABLED)
	if (proc_id == 0) {
                fclose (fout); 
                if (sim == 0)
                   fclose (fpvalue);
    	}
    	MPI_Finalize();
#	else
	fclose(fout);
	if(sim == 0) fclose(fpvalue);
#       endif

	currentCPUTime = clock();
	printf("time: %lf", (currentCPUTime - previousCPUTime) / (double) CLOCKS_PER_SEC);

  	return(1);
}

int Simulation (Tree *tree)
{
	int i, j, father;

	fprintf(fsim, "ID\tNUMBER");
	for(i=0; i<tree->ntaxa; i++)
		fprintf(fsim,"\t%s", tree->nodes[i].taxaname);
	fprintf(fsim, "\n");

	for(j=0; j<2*tree->ntaxa-1; j++){
		tree->nodes[j].contraction = 0;
		tree->nodes[j].nochange = 0;
		tree->nodes[j].expansion = 0;
	}

	for(i=0; i<ngenefamily; i++){
		tree->nodes[tree->root].theta=(int)((rndu()*(simngene[1]-simngene[0]+1))+simngene[0]);
		SimGeneFamily (tree->root, (int)tree->nodes[tree->root].theta, tree, simlambda[i]);
		//PrintTreeToFile(fout, &sptree);
		
		fprintf(fsim, "genefamily%d\tnumber%d", i+1, i+1);
		for(j=0; j<tree->ntaxa; j++) fprintf(fsim, "\t%d", (int)tree->nodes[j].theta);
		fprintf(fsim, "\n");

		for(j=0; j<2*tree->ntaxa-1; j++){
			if(j != tree->root){
				father = tree->nodes[j].father;
				if(tree->nodes[j].theta < tree->nodes[father].theta) tree->nodes[j].contraction++;
				else if (tree->nodes[j].theta == tree->nodes[father].theta) tree->nodes[j].nochange++;
				else tree->nodes[j].expansion++;
			}
		}
	}
	
	fprintf(fout, "lambda");
	for(i=0; i<2*tree->ntaxa-1; i++){
		if(i != tree->root) fprintf(fout,"contraction<node%d>\tnochange<node%d>\texpansion<node%d>\t",i+1,i+1,i+1);
	}
	fprintf(fout, "\n%lf\t",simlambda[0]);
	for(i=0; i<2*tree->ntaxa-1; i++){
		if(i != tree->root) fprintf(fout,"\t%d\t%d\t%d",tree->nodes[i].contraction, tree->nodes[i].nochange, tree->nodes[i].expansion);
	}
	fprintf(fout, "\n\n");

	return NO_ERROR;
}

int SimGeneFamily (int inode, int ncopies, Tree *tree, double lambda)
{
	int son1, son2;

	if(inode >= tree->ntaxa){
		son1 = tree->nodes[inode].sons[0];
		son2 = tree->nodes[inode].sons[1];
		tree->nodes[son1].theta = SimNcopies(ncopies, lambda, tree->nodes[son1].brlens);
		SimGeneFamily (son1, (int)tree->nodes[son1].theta, tree, lambda);
		tree->nodes[son2].theta = SimNcopies(ncopies, lambda, tree->nodes[son2].brlens);
 		SimGeneFamily (son2, (int)tree->nodes[son2].theta, tree, lambda);
	}
	return NO_ERROR;
}

int SimNcopies (int s, double lambda, double brlens)
{
	int ncopies, n;
	double prob, random;

	n = (int)(5*sqrt(2.0*s))+1;
	if(s < n) ncopies = 0;
	else ncopies = s - n;

	random = rndu();
	prob = LikelihoodBD (s, ncopies, lambda, brlens);
 
	while(random > prob){ 
		ncopies++;
		prob += LikelihoodBD(s, ncopies, lambda, brlens);
	}
            	
	return (ncopies);
}

void SimNcopiesProb (int s, double lambda, double brlens, int n, double *logprob)
{
	int i, j;
	int ncopies[100];
	double prob[100], total = 0.0, random;
	
	if(s < 50){
		for(i=0; i<100; i++) ncopies[i] = i;
	}else{
		for(i=0; i<100; i++) ncopies[i] = s - 50 + i + 1;
	}

	for(i=0; i<100; i++)	prob[i] = LikelihoodBD (s, ncopies[i], lambda, brlens);

	for(j=0; j<n; j++){
		random = rndu();
		for(i=0; i<100; i++){
			if(i == 0) total = prob[i];
			else total += prob[i];
			if(random <= total) break;
		}
		logprob[j] = log(prob[i]);
	}	 
}

int McmcRun (Tree *tree)
{
	int i, numnodes, numlambda;
	int round;
	int *numreject;
	double loglike, random, p;
	
	numreject = (int *)SafeMalloc((size_t) (ngenefamily * sizeof(int)));
	if (!numreject){
        	MrBayesPrint ("%s   Problem allocating numreject (%d)\n", spacer, ngenefamily * sizeof(int));
        	return (ERROR);
	}

	InitialParam (&sptree);

	/*initialize loglikelihood*/
	curLn = 0.0;
	for(i=0; i<ngenefamily; i++){
		if(Loglike1tree (tree, i, &loglike) == ERROR){
			printf("ERROR in LOGLIKE1TREE\n");
			return(ERROR);
		}
		curLn += loglike;
	}

	/*pick a move*/
	numnodes = ngenefamily * (tree->ntaxa - 1);
	numlambda = ngenefamily;

	for(round=1; round<mcmc.numGen; round++){
		random = rndu();
        	p = (double)numnodes/(numnodes+numlambda);

		if(random < p){
			MoveNode(tree);
		}else{
			MoveLambda(tree);
		}

		if(round % mcmc.sampleFreq == 0 || round == 1 || round == mcmc.numGen){
			PrintState (round, fout);
			if(Reject (numreject) == ERROR){
                      		printf("Error in Reject function\n");
                      		return(ERROR);
			}
			for(i=0; i<ngenefamily; i++) fprintf(fpvalue, "%d\t", numreject[i]);
            		fprintf(fpvalue, "\n");              
		}
	}	
	free (numreject);
	return (NO_ERROR);
}

int Reject (int *numreject)
{
    int i;
    double obsloglike, simloglike;

    for(i=0; i<ngenefamily; i++){
	     sptree.nodes[sptree.root].theta = sptree.nodes[sptree.root].ngenes[i];
	     SimGeneFamily (sptree.root, (int)sptree.nodes[sptree.root].theta, &sptree, sptree.nodes[0].lambda);

	     if(LoglikeSimtree (&sptree, &simloglike) == ERROR){
			printf("ERROR IN REJECTION");
			return (ERROR);
	     }

	     if(Loglike1tree (&sptree, i, &obsloglike) == ERROR){
			printf("ERROR IN REJECTION");
			return (ERROR);
	     }

             if(simloglike > obsloglike) numreject[i] = 0;
             else numreject[i] = 1;
    }

    return(NO_ERROR);   
}
                           
void InitialParam (Tree *tree)
{
	int i, j, k, n, sum_ngene;
	double maxbrlens=0.0, sum=0.0, sumsquare=0.0,t, lambda,x, variance;
	int *offspring;

	offspring = (int *)SafeMalloc((size_t) (tree->ntaxa * sizeof(int)));
	
	//initialize ngenes
	for(j=0; j<ngenefamily; j++){
		for(i=tree->ntaxa; i<2*(tree->ntaxa)-1; i++){
			//find descendant tips
			for(k=0; k<tree->ntaxa; k++) offspring[k] = 0;
			findOffsprings (offspring, tree, i);

			//mean number of gene copies
			sum_ngene = n = 0;
			for(k=0; k<tree->ntaxa; k++){
				if(offspring[k] == 1){
					sum_ngene += tree->nodes[k].ngenes[j];
					n++;
				}
			}
			tree->nodes[i].ngenes[j] = sum_ngene/n;
		}
	}

	//initialize lambda
	if(lambdalink == 1){
		t = TreeHeight(tree); 
		variance = 0.0;

		for(j=0; j<tree->ntaxa; j++){
			sum = sumsquare = 0.0;
			for(i=0; i<ngenefamily; i++){
				if(tree->nodes[tree->root].ngenes[i]==0) continue;
				x = tree->nodes[j].ngenes[i]/sqrt(2*(tree->nodes[tree->root].ngenes[i])*t); 
				sum += x;
				sumsquare += (x*x);
			}
			variance += (sumsquare - sum*sum/ngenefamily)/(ngenefamily-1);
		}

		lambda = variance/tree->ntaxa; 
		
		for(i=0; i<2*tree->ntaxa-1; i++){
			if(sptree.nodes[i].brlens > maxbrlens)
				maxbrlens = sptree.nodes[i].brlens;
		}

		if(lambda > 1/maxbrlens) lambda = 0.9/maxbrlens;

		for(i=0; i<2*tree->ntaxa-1; i++){
			tree->nodes[i].lambda = lambda;
			tree->nodes[i].maxlambda = 1.0/maxbrlens;
			tree->nodes[i].minlambda = 0.0;
			tree->nodes[i].lambdawindow = tree->nodes[i].maxlambda/5;
		}
	}else{
		t = TreeHeight(tree);
        
		for(i=0; i<2*tree->ntaxa-1; i++){
                        if(sptree.nodes[i].brlens > maxbrlens)
                        maxbrlens = sptree.nodes[i].brlens;                             
		}

                for(j=0; j<2*tree->ntaxa-1; j++){
                        sum = 0.0;
                        sumsquare = 0.0;
			tree->nodes[j].lambda = 0.0;

                        for(i=0; i<ngenefamily; i++){
                                if(tree->nodes[tree->root].ngenes[i]>0){
 					x = tree->nodes[j].ngenes[i]/sqrt(2*(tree->nodes[tree->root].ngenes[i])*t);
				}else{
					x = tree->nodes[j].ngenes[i]/sqrt(2*t);
				}
                                sum += x;
                                sumsquare += (x*x);
                        }

                        tree->nodes[j].lambda += (sumsquare - sum*sum/ngenefamily)/(ngenefamily-1); 
			if(tree->nodes[j].lambda > 1/maxbrlens) tree->nodes[j].lambda = 0.9/maxbrlens;
			tree->nodes[j].maxlambda = 1.0/maxbrlens;
                        tree->nodes[j].minlambda = 0.0;
                        tree->nodes[j].lambdawindow = tree->nodes[j].maxlambda/5;
                }

	}	
}

int ReadData (FILE *fin)
{
	int i, j, index = 0, *speciesindex;
	time_t t;
	struct tm *current;
	char datafile[30], outfile[30], simfile[30], pvaluefile[30];
	FILE *fdata; 
	char string[100], skip[100];

	fscanf(fin,"%d%s%ld%ld%d", &sim, datafile, &seed, &ngenefamily, &(sptree.ntaxa));


	if(sim == 0){
           sprintf(outfile, "%s.out", datafile);
           sprintf(pvaluefile,"%s.pvalue", datafile);
    	}else{
		sprintf(outfile, "%s.true", datafile);
		sprintf(simfile, "%s", datafile);
		fsim = (FILE*)gfopen(simfile,"w");
	}

#	if defined (MPI_ENABLED)
		if (proc_id == 0) 
        {
                    fout = (FILE*)gfopen(outfile,"w");
                    if(sim == 0) fpvalue = (FILE*)gfopen(pvaluefile,"w");
        }
#	else
		fout = (FILE*)gfopen(outfile,"w");
		if(sim == 0) fpvalue = (FILE*)gfopen(pvaluefile,"w");
#	endif

	/*set seed*/
	if(seed <= 0)
	{
		time(&t);
		current = localtime(&t);
#	if defined (MPI_ENABLED)
		seed = 11*current->tm_hour + 111*current->tm_min + (proc_id+1)*1111*current->tm_sec + 123;
#	else
		seed = 11*current->tm_hour + 111*current->tm_min + 1111*current->tm_sec + 123;
#	endif
		SetSeed(seed);
	}
	else
		SetSeed(seed);

	/*read the species tree*/
	ReadaTree(fin, &sptree);
	for(i=0; i<2*sptree.ntaxa-1; i++){
		sptree.nodes[i].ngenes = (int *)SafeMalloc((size_t) (ngenefamily * sizeof(int)));
		if (!sptree.nodes[i].ngenes){
			MrBayesPrint ("%s   Problem allocating sptree.nodes[%d].ngenes (%d)\n", spacer, i,ngenefamily * sizeof(int));
			return (ERROR);
		}
	}
	//PrintTreeToFile(fout, &sptree);

	/*simulation or MCMC*/
	if(sim == 1){
		simlambda = (double *)SafeMalloc((size_t)(ngenefamily * sizeof(double)));
		fscanf(fin,"%d%d", &simngene[0], &simngene[1]);
		for(i=0; i<ngenefamily; i++)
			fscanf(fin,"%lf", &simlambda[i]);
		return NO_ERROR;
	}else{
		fscanf(fin,"%d%d%d", &mcmc.numGen, &mcmc.sampleFreq, &lambdalink);
	}

	/*read gene family data*/
	fdata = (FILE*)gfopen(datafile,"r");
	speciesindex = (int *)SafeMalloc((size_t) (sptree.ntaxa * sizeof(int)));
	if (!speciesindex){
		MrBayesPrint ("Problem allocating speciesindex (%d)\n", ngenefamily * sizeof(int));
		return (ERROR);
	}
	for(i=0; i<sptree.ntaxa; i++) speciesindex[i] = 0;

	fscanf(fdata, "%s%s", skip, skip);
	for(i=0; i<sptree.ntaxa; i++){
		fscanf(fdata,"%s",string);
		for(j=0; j<sptree.ntaxa; j++){
			if(!strcmp(string,sptree.nodes[j].taxaname)){
				speciesindex[i] = j;
				break;
			}
		}
		if(speciesindex == 0){
			MrBayesPrint ("CANNOT FIND SPECIES %s in the species tree",string);
			return (ERROR);
		}
	}

	index = 0;
	while(!feof(fdata) & (index < ngenefamily)){
		fscanf(fdata, "%s%s",skip, skip);
		for(j=0; j<sptree.ntaxa; j++){
			fscanf(fdata, "%d", &(sptree.nodes[speciesindex[j]].ngenes[index]));
		}
		index++;
	}

	free (speciesindex);
	fclose(fdata);
	return (NO_ERROR);
}

void PrintHeader (void)
{
#		if !defined (MPI_ENABLED)
		MrBayesPrint ("\n\n");
#		endif
		MrBayesPrint ("%s            Bayesian estimation of Gene Duplication  \n\n",spacer);
#		if defined (MPI_ENABLED)
		MrBayesPrint ("                             (Parallel version)\n");
		MrBayesPrint ("                         (%d processors available)\n\n", num_procs);
#		endif
		srand ((unsigned int)time(NULL));
			MrBayesPrint ("%s                            by\n\n",spacer);
			MrBayesPrint ("%s                        Liang Liu\n\n",spacer);
			MrBayesPrint ("%s       Department of Agriculture and Natural Resources\n",spacer);
			MrBayesPrint ("%s                    Delaware State University\n",spacer);
			MrBayesPrint ("%s                         lliu@desu.edu\n\n",spacer);
		        MrBayesPrint ("%s       Distributed under the GNU General Public License\n\n",spacer);	
}




int PrintState (int round, FILE *outfile)
{
	char buffer[30];
  	struct timeval tv;
  	time_t curtime;
	int i, j, father,contraction, nochange, expansion;
	
	/*print to screen*/
	printf("%s round %d ---loglike: %f\n", spacer, round, curLn);
	
	/*print to file*/
	if(round == 1){
		gettimeofday(&tv, NULL);
        	curtime = tv.tv_sec;
               	strftime(buffer,30,"%T on %m-%d-%Y",localtime(&curtime));
		fprintf(outfile, "[This analysis was conducted at local time %s with tree ", buffer);
		PrintNodeToFile (outfile, &sptree);
		fprintf(outfile, ". The numbers in the tree are node numbers.]\n");

		fprintf(outfile,"loglike\t");
		if(lambdalink == YES){
			fprintf(outfile,"lambda\t");
		}else{
			for(i=0; i<2*sptree.ntaxa-1; i++)
			{
			if(i == sptree.root) continue;
			fprintf(outfile,"lambda<node%d>\t", i+1);
			}
		}

		for(i=0; i<2*sptree.ntaxa-1; i++){
			if(i != sptree.root)
				fprintf(outfile,"contraction<node%d>\tnochange<node%d>\texpansion<node%d>\t",i+1,i+1,i+1);
		}

		fprintf(outfile,"\n");
		fflush(outfile);
	}else{
		fprintf(outfile,"%lf\t",curLn);
		if(lambdalink == YES){
			fprintf(outfile,"%lf\t",sptree.nodes[0].lambda);
		}else{
			for(i=0; i<2*sptree.ntaxa-1; i++){
				if(i == sptree.root) continue;
				fprintf(outfile,"%lf\t",sptree.nodes[i].lambda);
			}
		}

		for(i=0; i<2*sptree.ntaxa-1; i++){
			if(i != sptree.root){
				father = sptree.nodes[i].father;
				contraction = 0;
				expansion = 0;
				nochange = 0;
				for(j=0; j<ngenefamily; j++){
					if(sptree.nodes[i].ngenes[j] < sptree.nodes[father].ngenes[j]) contraction++;
					else if (sptree.nodes[i].ngenes[j] == sptree.nodes[father].ngenes[j]) nochange++;
					else expansion++;
				}
				fprintf(outfile,"%d\t%d\t%d\t",contraction, nochange, expansion);
			}
		}
		
		fprintf(outfile,"\n");
		fflush(outfile);
	}

	return (NO_ERROR);
}

int MoveNode (Tree *tree)
{
	int treeindex, nodeindex, oldpara, newpara;
	double x, oldloglike, newloglike, diff, random;

	/*pick a node at random*/
	treeindex = rndu() * ngenefamily;
	nodeindex = rndu() * (tree->ntaxa - 1) + tree->ntaxa;
	oldpara = tree->nodes[nodeindex].ngenes[treeindex];

	/*calculate likelihood*/
	if(Loglike1tree (tree, treeindex, &oldloglike) == ERROR){
		printf("Error in Loglike1tree\n");
		return (ERROR);
	}

	/*new parameter value*/
	if(oldpara == 0){
		x = rndu();
		if(x < 0.5)
			newpara = oldpara;
		else
			newpara = oldpara + 1;
	}else{
		x = rndu();
		if (x < 0.333333333)
			newpara = oldpara + 1;
		else if (0.333333333 < x && x < 0.666666666)
			newpara = oldpara - 1;
		else
			newpara = oldpara;
	}
	tree->nodes[nodeindex].ngenes[treeindex] = newpara;

	/*new loglikelihood*/	
	if(Loglike1tree (tree, treeindex, &newloglike) == ERROR){
		printf("Error in Loglike1tree\n");
		return (ERROR);
	}
	
	/*move*/
	diff = newloglike - oldloglike;
	random = log(rndu ());
	if(random > diff)
		tree->nodes[nodeindex].ngenes[treeindex] = oldpara;
	else
		curLn += diff;

	return (NO_ERROR);
}

int MoveLambda (Tree *tree)
{
	int i, nodeindex;
	double max, min, window, oldlambda, newlambda, loglike, oldloglike, newloglike, diff, random;

	if(lambdalink == YES){
		window = tree->nodes[0].lambdawindow;
		oldlambda = tree->nodes[0].lambda;
		oldloglike = curLn;

		max = Min(tree->nodes[0].maxlambda, oldlambda+window);
		min = Max(tree->nodes[0].minlambda, oldlambda-window);
	
		newlambda = rndu() * (max-min) + min;
		for(i=0; i<2*tree->ntaxa-1; i++)
			tree->nodes[i].lambda = newlambda;


		newloglike = 0.0;
		for(i=0; i<ngenefamily; i++){
			if(Loglike1tree (tree, i, &loglike) == ERROR){
				printf("ERROR in LOGLIKE1TREE\n");
				return(ERROR);
			}
			newloglike += loglike;
		}

		diff = newloglike - oldloglike;
		random = log(rndu ());
		if(random > diff)
			for(i=0; i<2*tree->ntaxa-1; i++) tree->nodes[i].lambda = oldlambda;
		else curLn += diff;
	
	}else{
		/*pick a node at random*/
		do{
			nodeindex = rndu() * (2*tree->ntaxa - 1);
		}while(nodeindex == tree->root);

		window = tree->nodes[nodeindex].lambdawindow;
		oldlambda = tree->nodes[nodeindex].lambda;
		oldloglike = curLn;

		max = Min(tree->nodes[nodeindex].maxlambda, oldlambda+window);
		min = Max(tree->nodes[nodeindex].minlambda, oldlambda-window);
		
		newlambda = rndu() * (max-min) + min;
		tree->nodes[nodeindex].lambda = newlambda;

		newloglike = 0.0;
		for(i=0; i<ngenefamily; i++){
			if(Loglike1tree (tree, i, &loglike) == ERROR){
				printf("ERROR in LOGLIKE1TREE\n");
				return(ERROR);
			}
			newloglike += loglike;
		}

		diff = newloglike - oldloglike;
		random = log(rndu ());
		if(random > diff) tree->nodes[nodeindex].lambda = oldlambda;
		else curLn += diff;
	}

	return (NO_ERROR);
}

double LikelihoodBD (int s, int c, double lambda, double brlens)
{
	int j;
	int m = Min(s,c);
	double alpha, likelihood = 0.0;

	
	if(lambda > 1.0/brlens || (s == 0 && c > 0)){
		return (0.0);
	}else if(s == 0 && c == 0){
		return (1.0);
	}else{
		alpha = lambda*brlens/(1+lambda*brlens);

		for(j=0; j<=m; j++){
			likelihood += exp(lnchoose(s,j) + lnchoose(s+c-j-1,s-1) + log(alpha) * (s+c-2*j) + log(1-2*alpha)*j);
		}
		return (likelihood);
	}
}

int Loglike1tree (Tree *tree, int igene, double *loglike)
{
	int i, s, c, father;
	double brlens, lambda, likelihood;

	*loglike = 0.0;
	for(i=0; i<2*(tree->ntaxa)-1; i++){
		if(i == tree->root){
 			continue;
		}else{
			father = tree->nodes[i].father; 
			s = tree->nodes[father].ngenes[igene];
			c = tree->nodes[i].ngenes[igene];
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			
			likelihood = LikelihoodBD (s, c, lambda, brlens);

			if(likelihood > 0){
 				*loglike += log(likelihood);
			}else{
				*loglike = -99999.9;
				break;
			}
		}
	}
	return (NO_ERROR);
}


int LoglikeSimtree (Tree *tree, double *loglike)
{
	int i, s, c, father;
	double brlens, lambda, likelihood;

	
	*loglike = 0.0;
	for(i=0; i<2*(tree->ntaxa)-1; i++){
		if(i == tree->root) continue;
		else{
			father = tree->nodes[i].father; 
			s = tree->nodes[father].theta;
			c = (int)tree->nodes[i].theta;
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			
			likelihood = LikelihoodBD (s, c, lambda, brlens);

			if(likelihood > 0) *loglike += log(likelihood);
			else
			{
				*loglike = -99999.9;
				break;
			}
		}
	}
	return (NO_ERROR);
}











