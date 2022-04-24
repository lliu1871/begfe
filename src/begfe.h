#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>
#include <stdarg.h>
#include <sys/time.h>


#define Max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#define Min( a, b ) ( ((a) < (b)) ? (a) : (b) )

#define NTAXA         100      /* max # of species */
#define NGENE         50000      /* max # of loci */
#define MAXROUND	1000000		/* MAX # OF ROUNDS*/
#define NUM_NOCHANGE	10000		/* # OF ROUNDS THAT NO BIGGER LIKELIHOOD VALUES ARE FOUND*/
#define LSPNAME       30       /* # characters in sequence names */
#define ERROR 1
#define NO_ERROR 0
#define YES 1
#define NO 0
#define NA -1
#define DEBUG 0
#define FPN(file) fputc('\n', file)
#define FOR(i,n) for(i=0; i<n; i++)
#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))
#define CDFGamma(x,alpha,beta) IncompleteGamma((beta)*(x),alpha,LnGamma(alpha))
typedef struct node 
	{
	int		father, nson, sons[2], namenumber;
	char	taxaname[LSPNAME];
	int		*ngenes, expansion, contraction, nochange;
	double	brlens, theta, lambda;
   	}
	Treenode;
typedef struct Tree
	{
   	int root;
	int ntaxa; 
   	Treenode nodes[2*NTAXA];
	}  
	Tree;
typedef struct chain
	{
	int			numGen;                /* number of MCMC cycles                         */
	int			sampleFreq;            /* frequency to sample chain                     */
	int			printFreq;             /* frequency to print chain                      */
	int			swapFreq;              /* frequency to attempt swap of states           */
	int			numRuns;               /* number of runs                                */
	int			numChains;             /* number of chains                              */
	} Chain;
