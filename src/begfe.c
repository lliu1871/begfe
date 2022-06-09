/*
 *  Begfe version 2.0
 *
 *  copyright 2022-2025
 *
 *  Liang Liu
 *  Department of Statistics and Institute of Bioinformatics
 *  University of Georgia
 *
 *  lliu@uga.edu
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

#if defined(PTHREAD_ENABLED)
#include 	<pthread.h>
#define 	NUMTHREADS 8
#endif


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
int 		SimGeneFamily (int inode, int ncopies, Tree *tree);
int 		SimNcopies (int s, double lambda, double brlens);
int			Simulation (Tree *tree);
int 		Loglike1treeThread (Tree *tree, int igene, double *loglike);

/*local variables*/
Tree 		sptree;
double		maxlambda=1.0, minlambda=0.0, lambda_window=0.1;
long int	seed=0;
long int	ngenefamily;
double		curLn, *curLnGenetrees;
FILE		*fout;
FILE        *fpvalue;
Chain	 	mcmc;
int			lambdalink;
int			sim;
int			simngene[2];

int main (int argc, char *argv[])
{
	FILE *fin;

	PrintHeader();

	fin = (FILE*)gfopen(argv[1],"r");
	if(ReadData (fin) == ERROR){
		MrBayesPrint ("Problem in ReadData\n");
		return ERROR;
	}
	fclose(fin);	

	if(sim == 1) {
		if(Simulation (&sptree) == ERROR){
			MrBayesPrint ("Errors in MCMCRUN\n");
			return ERROR;
		}
	}else{
		if(McmcRun (&sptree) == ERROR){
			MrBayesPrint ("Errors in MCMCRUN\n");
			return ERROR;
		}
	}
	
	/*free memory*/
	if(sim == 0){ 
		fclose(fpvalue);
		fclose(fout);
	}else{
		free(curLnGenetrees);
		fclose(fout);
	}
  	return(NO_ERROR);
}

int Simulation (Tree *tree)
{
	int i, j;

	fprintf(fout, "ID");
	for(i=0; i<tree->ntaxa; i++){
		fprintf(fout,"\t%s", tree->nodes[i].taxaname);
	}

	for(j=0; j<2*tree->ntaxa-1; j++){
		tree->nodes[j].contraction = 0;
		tree->nodes[j].nochange = 0;
		tree->nodes[j].expansion = 0;
	}

	for(i=0; i<ngenefamily; i++){
		tree->nodes[tree->root].theta=(int)((rndu()*(simngene[1]-simngene[0]+1))+simngene[0]);
		SimGeneFamily (tree->root, (int)tree->nodes[tree->root].theta, tree);
		//PrintTreeToFile(fout, &sptree);
		
		fprintf(fout, "\ngenefamily%d\t", i+1);
		for(j=0; j<tree->ntaxa; j++){ 
			fprintf(fout, "\t%d", (int)tree->nodes[j].theta);
		}
	}
	return NO_ERROR;
}

int SimGeneFamily (int inode, int ncopies, Tree *tree)
{
	int son1, son2;

	if(inode >= tree->ntaxa){
		son1 = tree->nodes[inode].sons[0];
		son2 = tree->nodes[inode].sons[1];
		tree->nodes[son1].theta = SimNcopies(ncopies, tree->nodes[inode].lambda, tree->nodes[son1].brlens);
		SimGeneFamily (son1, (int)tree->nodes[son1].theta, tree);
		tree->nodes[son2].theta = SimNcopies(ncopies, tree->nodes[inode].lambda, tree->nodes[son2].brlens);
 		SimGeneFamily (son2, (int)tree->nodes[son2].theta, tree);
	}
	
	return NO_ERROR;
}

int SimNcopies (int s, double lambda, double brlens)
{
	int ncopies = 0;
	double prob, random;
	
	/*if s = 0, then ncopies = 0*/
	if(s == 0){
		return(ncopies);
	}

	/*calculate the cumulative birth-death probability*/
	random = rndu();
	prob = LikelihoodBD (s, ncopies, lambda, brlens);	
	while(random > prob){ 
		ncopies++;
		prob += LikelihoodBD(s, ncopies, lambda, brlens);
		
		/*to avoid dead loop*/
		if(ncopies > s*100){
			printf("The total probability is not 1.0 when ncopies = %d\n", ncopies);
			break;
		}
	}	
	return (ncopies);
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
		curLnGenetrees[i] = loglike;
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
			/*print out posterior in file*/
			PrintState (round, fout);
			
			/*calculate Bayesian posterior pvalue*/
			if(Reject (numreject) == ERROR){
				printf("Error in Reject function\n");
				return(ERROR);
			}		
			for(i=0; i<ngenefamily; i++){ 
				fprintf(fpvalue, "%d\t", numreject[i]);
			}
            fprintf(fpvalue, "\n");   
			fflush(fpvalue);
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
		 SimGeneFamily (sptree.root, (int)sptree.nodes[sptree.root].theta, &sptree);
	     if(LoglikeSimtree (&sptree, &simloglike) == ERROR){
			printf("ERROR IN REJECTION");
			return (ERROR);
	     }

	     if(Loglike1tree (&sptree, i, &obsloglike) == ERROR){
			printf("ERROR IN REJECTION");
			return (ERROR);
	     }

         if(simloglike > obsloglike){ 
			numreject[i] = 0;
		 }else{ 
			numreject[i] = 1;
		 }
    }

    return(NO_ERROR);   
}
                           
void InitialParam (Tree *tree)
{
	int i, j, k, n, sum_ngene;
	double sum=0.0, sumsquare=0.0, t, lambda, x, variance;
	int *offspring;

	//initialize curLnGenetrees
	curLnGenetrees	= (double *)SafeMalloc((size_t) (ngenefamily * sizeof(double)));
	
	//initialize # of gene copies
	offspring = (int *)SafeMalloc((size_t) (tree->ntaxa * sizeof(int)));
	for(j=0; j<ngenefamily; j++){
		for(i=tree->ntaxa; i<2*(tree->ntaxa)-1; i++){
			//find descendant tips
			for(k=0; k<tree->ntaxa; k++){ 
				offspring[k] = 0;
			}
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
			if(tree->nodes[i].ngenes[j] < 1) tree->nodes[i].ngenes[j] += 1;
		}
	}
	free(offspring);
	
	//initialize lambda
	t = TreeHeight(tree); 
	variance = 0.0;

	/*if ngenefamily > 1, lambda is equal to the average variance*/
	if(ngenefamily > 1){
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
	}

	/*if ngenefamily = 1, lambda is the variance of gene copy numbers*/	
	if(ngenefamily == 1){
		sum = sumsquare = 0.0;
		for(j=0; j<tree->ntaxa; j++){
			x = tree->nodes[j].ngenes[0]/sqrt(2*(tree->nodes[tree->root].ngenes[0])*t); 
			sum += x;
			sumsquare += (x*x);
		}
		variance = sumsquare/(tree->ntaxa-1) - sum*sum/(tree->ntaxa*(tree->ntaxa-1));
		lambda = variance;
	}

	if(lambda > maxlambda){
		lambda = 0.9 * maxlambda;
	}

	for(i=0; i<2*tree->ntaxa-1; i++){
		tree->nodes[i].lambda = lambda;
	}
}

int ReadData (FILE *fin)
{
	int i, j, index = 0, *speciesindex;
	time_t t;
	struct tm *current;
	char datafile[30], outfile[50], simfile[50], pvaluefile[50];
	FILE *fdata; 
	char string[100], skip[100];

	if(fscanf(fin,"%d%s%ld%ld%d", &sim, datafile, &seed, &ngenefamily, &(sptree.ntaxa)) != 5){
		MrBayesPrint ("%s   Problem reading the control file\n");
		return ERROR;
	}

	if(sim == 0){
		sprintf(outfile, "%s.out", datafile);
		sprintf(pvaluefile,"%s.pvalue", datafile);
		fpvalue = (FILE*)gfopen(pvaluefile,"w");
		fout = (FILE*)gfopen(outfile,"w");
    }else{
		sprintf(simfile, "%s", datafile);
		fout = (FILE*)gfopen(simfile,"w");
	}

	/*set seed*/
	if(seed <= 0){
		time(&t);
		current = localtime(&t);
		seed = 11*current->tm_hour + 111*current->tm_min + 1111*current->tm_sec + 123;
		SetSeed(seed);
	}else{
		SetSeed(seed);
	}

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
		if(fscanf(fin,"%d%d", &simngene[0], &simngene[1]) != 2){
			MrBayesPrint ("%s   Problem reading the control file\n");
			return ERROR;
		}else{
			return NO_ERROR;
		}
	}else{
		if(fscanf(fin,"%d%d%d", &mcmc.numGen, &mcmc.sampleFreq, &lambdalink) != 3){
			MrBayesPrint ("%s   Problem reading the control file\n");
			return ERROR;
		}
	}

	/*read gene family data*/
	fdata = (FILE*)gfopen(datafile,"r");
	speciesindex = (int *)SafeMalloc((size_t) (sptree.ntaxa * sizeof(int)));
	if (!speciesindex){
		MrBayesPrint ("Problem allocating speciesindex (%d)\n", ngenefamily * sizeof(int));
		return (ERROR);
	}
	for(i=0; i<sptree.ntaxa; i++) speciesindex[i] = -1;

	if(!fscanf(fdata, "%s", skip)){
		MrBayesPrint ("%s   Problem reading the control file\n");
		return ERROR;
	}
		
	for(i=0; i<sptree.ntaxa; i++){
		if(!fscanf(fdata,"%s",string)){
			MrBayesPrint ("%s   Problem reading the control file\n");
			return ERROR;
		}
		for(j=0; j<sptree.ntaxa; j++){
			if(!strcmp(string,sptree.nodes[j].taxaname)){
				speciesindex[i] = j;
				break;
			}
		}
		if(speciesindex[i] == -1){
			MrBayesPrint ("Cannot find species %s in the species tree\n",string);
			return (ERROR);
		}
	}

	index = 0;
	while(!feof(fdata) & (index < ngenefamily)){
		if(!fscanf(fdata, "%s",skip)){
			MrBayesPrint ("%s   Problem reading the control file\n");
			return ERROR;
		}
		for(j=0; j<sptree.ntaxa; j++){
			if(!fscanf(fdata, "%d", &(sptree.nodes[speciesindex[j]].ngenes[index]))){
				MrBayesPrint ("%s   Problem reading the control file\n");
				return ERROR;
			}
		} 
		index++;
	}

	if(index < ngenefamily){
		printf("# of gene families %d != %ld\n", index, ngenefamily);
		return(ERROR);
	}
	
	free (speciesindex);
	fclose(fdata);
	return (NO_ERROR);
}

void PrintHeader (void)
{
	MrBayesPrint ("%s       Bayesian Estimation of Gene Duplication (version 2.0)  \n\n",spacer);
	MrBayesPrint ("%s                        Liang Liu\n",spacer);
	MrBayesPrint ("%s                   University of Georgia\n",spacer);
	MrBayesPrint ("%s                       lliu@uga.edu\n\n",spacer);
	MrBayesPrint ("%s         Distributed under the GNU General Public License\n\n",spacer);	
}

int PrintState (int round, FILE *outfile)
{
	char buffer[30];
  	struct timeval tv;
  	time_t curtime;
	int i, j, father,contraction, nochange, expansion;
	
	/*print to screen*/
	MrBayesPrint ("%s         round %d --------------- loglike: %f\n", spacer, round, curLn);
	
	/*print to file*/
	if(round == 1){
		gettimeofday(&tv, NULL);
        curtime = tv.tv_sec;
        strftime(buffer,30,"%X on %m-%d-%Y",localtime(&curtime));
		fprintf(outfile, "[This analysis was conducted at local time %s with seed = %ld and tree ", buffer, seed);
		PrintNodeToFile (outfile, &sptree);
		fprintf(outfile, ". The numbers in the tree are node numbers.]\n");

		fprintf(outfile,"round\tloglike\t");
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
	}
	
	fprintf(outfile,"\n%d\t%lf\t",round, curLn);
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
	fflush(outfile);
	
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
	oldloglike = curLnGenetrees[treeindex];

	/*new parameter value*/
	if(oldpara == 1){
		x = rndu();
		if(x < 0.5)
			newpara = oldpara;
		else
			newpara = oldpara + 1;
	}else{
		x = rndu();
		if (x < 0.333333333){
			newpara = oldpara + 1;
		}else if (0.333333333 < x && x < 0.666666666){
			newpara = oldpara - 1;
		}else{
			newpara = oldpara;
		}
	}
	
	tree->nodes[nodeindex].ngenes[treeindex] = newpara;

	/*new loglikelihood*/	
	if(Loglike1tree (tree, treeindex, &newloglike) == ERROR){
		printf("Error in Loglike1tree\n");
		return (ERROR);
	}
	
	/*move*/
	diff = newloglike - oldloglike;
	random = log(rndu());
	if(random > diff){
		tree->nodes[nodeindex].ngenes[treeindex] = oldpara;
	}else{
		curLn += diff;
		curLnGenetrees[treeindex] = newloglike;
	}

	return (NO_ERROR);
}

int MoveLambda (Tree *tree)
{
	int i, nodeindex;
	double max, min, oldlambda, newlambda, loglike, oldloglike, newloglike, newloglikegenetree[NGENEFAMILY], diff, random;

	if(lambdalink == YES){
		oldlambda = tree->nodes[0].lambda;
		oldloglike = curLn;

		max = Min(maxlambda, oldlambda + lambda_window);
		min = Max(minlambda, oldlambda - lambda_window);
	
		newlambda = rndu() * (max-min) + min;
		for(i=0; i<2*tree->ntaxa-1; i++){
			tree->nodes[i].lambda = newlambda;
		}

		newloglike = 0.0;
		for(i=0; i<ngenefamily; i++){
			if(Loglike1tree (tree, i, &loglike) == ERROR){
				printf("ERROR in LOGLIKE1TREE\n");
				return(ERROR);
			}
			newloglikegenetree[i] = loglike;
			newloglike += loglike;
		}

		diff = newloglike - oldloglike;
		random = log(rndu ());
		if(random > diff){
			for(i=0; i<2*tree->ntaxa-1; i++){
				tree->nodes[i].lambda = oldlambda;
			}
		}else{ 
			/*update curLn and curLnGenetrees*/
			curLn += diff;
			for(i=0; i<ngenefamily; i++){
				curLnGenetrees[i] = newloglikegenetree[i];
			}
		}
	
	}else{
		/*pick a node (excluding the tree root) at random*/
		do{
			nodeindex = rndu() * (2*tree->ntaxa - 1);
		}while(nodeindex == tree->root);

		oldlambda = tree->nodes[nodeindex].lambda;
		oldloglike = curLn;

		max = Min(maxlambda, oldlambda + lambda_window);
		min = Max(minlambda, oldlambda - lambda_window);
		
		newlambda = rndu() * (max-min) + min;
		tree->nodes[nodeindex].lambda = newlambda;

		newloglike = 0.0;
		for(i=0; i<ngenefamily; i++){
			if(Loglike1tree (tree, i, &loglike) == ERROR){
				printf("ERROR in LOGLIKE1TREE\n");
				return(ERROR);
			}
			newloglikegenetree[i] = loglike;
			newloglike += loglike;
		}

		diff = newloglike - oldloglike;
		random = log(rndu ());
		if(random > diff){ 
			tree->nodes[nodeindex].lambda = oldlambda;
		}else{ 
			curLn += diff;
			for(i=0; i<ngenefamily; i++){
				curLnGenetrees[i] = newloglikegenetree[i];
			}
		}
	}

	return (NO_ERROR);
}

double LikelihoodBD (int s, int c, double lambda, double brlens)
{
	int j;
	int m = Min(s,c);
	double alpha, likelihood;

	alpha = lambda*brlens/(1+lambda*brlens);
	
	if(s == 0 && c > 0){
		likelihood = 0.0;
		printf("Error: s=0 and c>0\n");
		exit(-1);
	}else if(s == 0 && c == 0){
		likelihood = 1.0;
	}else if(s > 0 && c == 0){
		likelihood = exp(s*log(alpha));
	}else{
		likelihood = 0.0;
		
		if(alpha < 0.5){
			for(j=0; j<=m; j++){
				likelihood += exp(lnchoose(s,j) + lnchoose(s+c-j-1,s-1) + log(alpha) * (s+c-2*j) + log(1-2*alpha)*j);
			}
		}else{
			for(j=1; j<=s; j++){
				likelihood += exp(lnchoose(s,j) + lnchoose(c+j-1,j-1) + log(2*alpha-1) * (s-j) + log(1-alpha)*2*j + log(alpha)*(c-s));
			}
		}
		
	}
	return (likelihood);
}


#ifdef NUMTHREADS
struct ThreadData {
    int start, stop, igene, threadID;
    Tree *tree;
	double *loglikelihood;
};

void *loglikeThread (void *td) {
    struct ThreadData* data=(struct ThreadData*) td;
    int start = data->start;
    int stop = data->stop;	
	int igene = data->igene;
	int threadID = data->threadID;	
    double* loglike = data->loglikelihood;
	Tree* tree = data->tree;	
	int i, father, s, c;
	double brlens, lambda, likelihood;

	loglike[threadID] = 0.0;
    for (i=start; i<stop; i++) {
        if(i != tree->root){
			father = tree->nodes[i].father; 
			s = tree->nodes[father].ngenes[igene];
			c = tree->nodes[i].ngenes[igene];
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			likelihood = LikelihoodBD (s, c, lambda, brlens);
 			loglike[threadID] += log(likelihood);
		}
    }
    return NULL;
}

int Loglike1tree (Tree *tree, int igene, double *loglike)
{
	int i;
	double loglikelihood[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    struct ThreadData data[NUMTHREADS];

	/*
        this has the effect of rounding up the number of tasks
        per thread, which is useful in case ARRAYSIZE does not
        divide evenly by NUMTHREADS.
    */
    int tasksPerThread=(2*(tree->ntaxa)-1)/NUMTHREADS;
	
    /* Divide work for threads, prepare parameters */
    for (i=0; i<NUMTHREADS; i++) {
        data[i].start=i*tasksPerThread;
        data[i].stop=(i+1)*tasksPerThread;
        data[i].tree=tree;
		data[i].loglikelihood=loglikelihood;
		data[i].threadID = i;
		data[i].igene=igene;
    }
	
    /* the last thread must not go past the end of the array */
    data[NUMTHREADS-1].stop = 2*(tree->ntaxa)-1;
	
    /* Launch Threads */
    for (i=0; i<NUMTHREADS; i++) {
        pthread_create(&thread[i], NULL, loglikeThread, (void*)&data[i]);
    }

    /* Wait for Threads to Finish */
    for (i=0; i<NUMTHREADS; i++) {
        pthread_join(thread[i], NULL);
    }

	*loglike = 0.0;
    /* Display Result */
    for (i=0; i<NUMTHREADS; i++) {
		*loglike += loglikelihood[i];
    }
	
	return (NO_ERROR);
}

void *loglikeSimThread (void *td) {
    struct ThreadData* data=(struct ThreadData*) td;
    int start = data->start;
    int stop = data->stop;	
	//int igene = data->igene;
	int threadID = data->threadID;	
    double* loglike = data->loglikelihood;
	Tree* tree = data->tree;	
	int i, father, s, c;
	double brlens, lambda, likelihood;

	loglike[threadID] = 0.0;
    for (i=start; i<stop; i++) {
        if(i != tree->root){
			father = tree->nodes[i].father; 
			s = tree->nodes[father].theta;
			c = tree->nodes[i].theta;
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			likelihood = LikelihoodBD (s, c, lambda, brlens);
 			loglike[threadID] += log(likelihood);
		}
    }
    return NULL;
}

int LoglikeSimtree (Tree *tree, double *loglike)
{
	int i;
	double loglikelihood[NUMTHREADS];
    pthread_t thread[NUMTHREADS];
    struct ThreadData data[NUMTHREADS];

	/*
        this has the effect of rounding up the number of tasks
        per thread, which is useful in case ARRAYSIZE does not
        divide evenly by NUMTHREADS.
    */
    int tasksPerThread=(2*(tree->ntaxa)-1)/NUMTHREADS;

    /* Divide work for threads, prepare parameters */
    for (i=0; i<NUMTHREADS; i++) {
        data[i].start=i*tasksPerThread;
        data[i].stop=(i+1)*tasksPerThread;
        data[i].tree=tree;
		data[i].loglikelihood=loglikelihood;
		data[i].threadID = i;
    }
	
    /* the last thread must not go past the end of the array */
    data[NUMTHREADS-1].stop = 2*(tree->ntaxa)-1;
	
    /* Launch Threads */
    for (i=0; i<NUMTHREADS; i++) {
        pthread_create(&thread[i], NULL, loglikeSimThread, (void*)&data[i]);
    }

    /* Wait for Threads to Finish */
    for (i=0; i<NUMTHREADS; i++) {
        pthread_join(thread[i], NULL);
    }

	*loglike = 0.0;
    /* Display Result */
    for (i=0; i<NUMTHREADS; i++) {
		*loglike += loglikelihood[i];
    }
	
	return (NO_ERROR);
}

#else

int Loglike1tree (Tree *tree, int igene, double *loglike)
{
	int i, s, c, father;
	double brlens, lambda, likelihood;

	*loglike = 0.0;
	for(i=0; i<2*(tree->ntaxa)-1; i++){
		if(i != tree->root){
			father = tree->nodes[i].father; 
			s = tree->nodes[father].ngenes[igene];
			c = tree->nodes[i].ngenes[igene];
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			likelihood = LikelihoodBD (s, c, lambda, brlens);
 			*loglike += log(likelihood);
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
		if(i != tree->root){
			father = tree->nodes[i].father; 
			s = tree->nodes[father].theta;
			c = (int)tree->nodes[i].theta;
			brlens = tree->nodes[i].brlens;
			lambda = tree->nodes[i].lambda;
			likelihood = LikelihoodBD (s, c, lambda, brlens);
			*loglike += log(likelihood);
		}
	}
	return (NO_ERROR);
}
#endif	






