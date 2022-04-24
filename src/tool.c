#include "begfe.h"
#include "tool.h"

int 		AddToPrintString (char *tempStr);
int 		SaveSprintf(char **target, int *targetLen, char *fmt, ...);
void		MrBayesPrint (char *format, ...);
int 		PrintState (int round, FILE *fout, int addend);
int 		PrintNode(Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);
int 		PrintTree (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void 		WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted);
void 		WriteNodeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted);

char		spacer[10]="  ";
char		*printString;                /* string for printing to a file                */
size_t		printStringSize;             /* length of printString                        */

FILE *gfopen(char *filename, char *mode)
{
   FILE *fp=(FILE*)fopen(filename, mode);
   if(fp==NULL) {
      printf("\nerror when opening file %s\n", filename);
      exit(-1);
   }
   return(fp);
}

void error2 (char * message)
{ printf("\nError: %s.\n", message); exit(-1); }


static time_t time_start;
void starttime (void)
{  
   time_start=time(NULL);
}




void strcase (char *str, int direction)
{
/* direction = 0: to lower; 1: to upper */
   char *p=str;
   if(direction)  while(*p) { *p=(char)toupper(*p); p++; }
   else           while(*p) { *p=(char)tolower(*p); p++; }
}


int matIout (FILE *fout, int x[], int n, int m)
{
   int i,j;
   for (i=0,FPN(fout); i<n; i++,FPN(fout)) 
      FOR(j,m) fprintf(fout,"%6d", x[i*m+j]);
   return (0);
}


/*///////////////////////////////////////////*/
/*      math function                       */
/*/////////////////////////////////////////*/
long factorial(int n)
{
   long f, i;
   if (n>10) error2("n>10 in factorial");
   for (i=2,f=1; i<=(long)n; i++) f*=i;
   return (f);
}

double lnchoose (int n, int x)
{
	int i;
	double f=0.0;

	if(n < x){
		printf("n must be greater than or equal to x\n");
		exit(-1);
	}else if (x == 0 || n == x){
		return 0.0;
	}else {
		for(i=n; i >= n-x+1; i--)
			f += (log(i)-log(n-i+1));
		return f;
	}
}

static unsigned int z_rndu=137;
static unsigned int w_rndu=123456757;

void SetSeed (unsigned int seed)
{
   if(sizeof(int)-4) puts("oh-oh, we are in trouble.  int not 32-bit?");
   z_rndu=170*(seed%178)+137;
   w_rndu = seed*127773;
}

double rndu (void)
{
/* U(0,1): AS 183: Appl. Stat. 31:188-190 
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190


   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
   static unsigned int x_rndu=11, y_rndu=23;
   double r;
   
   x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
   y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
   z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
/*
   if (x_rndu<0) x_rndu+=30269;
   if (y_rndu<0) y_rndu+=30307;
   if (z_rndu<0) z_rndu+=30323;
*/
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;


  return (r-(int)r);
}


double rndgamma1 (double s);
double rndgamma2 (double s);


double rndgamma (double s)
{
/* random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
      r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
           rndgamma2() if s>1 or
           exponential if s=1
*/
   double r=0;

   if (s<=0)      puts ("jgl gamma..");
   else if (s<1)  r=rndgamma1 (s);
   else if (s>1)  r=rndgamma2 (s);
   else           r=-log(rndu());
   return (r);
}




double rndgamma1 (double s)
{
/* random standard gamma for s<1
   switching method
*/
   double r, x=0,small=1e-37,w;
   static double a,p,uf,ss=10,d;


   if (s!=ss) {
      a=1-s;
      p=a/(a+s*exp(-a));
      uf=p*pow(small/a,s);
      d=a*log(a);
      ss=s;
   }
   for (;;) {
      r=rndu();
      if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
      else if (r>uf)  x=a*pow(r/p,1/s), w=x;
      else            return (0);
      r=rndu ();
      if (1-r<=w && r>0)
         if (r*(w+1)>=1 || -log(r)<=w)  continue;
      break;
   }
   return (x);
}


double rndgamma2 (double s)
{
/* random standard gamma for s>1
   Best's (1978) t distribution method
*/
   double r,d,f,g,x;
   static double b,h,ss=0;
   if (s!=ss) {
      b=s-1;
      h=sqrt(3*s-0.75);
      ss=s;
   }
   for (;;) {
      r=rndu ();
      g=r-r*r;
      f=(r-0.5)*h/sqrt(g);
      x=b+f;
      if (x <= 0) continue;
      r=rndu();
      d=64*r*r*g*g*g;
      if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
   }
   return (x);
}


int comparedouble (const void *a, const void *b)
{  
   double aa = *(double*)a, bb= *(double*)b;
   return (aa==bb ? 0 : aa > bb ? 1 : -1);
}


double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
           limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   /* double accurate=1e-8, overflow=1e30; */
   double accurate=1e-10, overflow=1e60;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];


   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);


   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;


   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;


 l50:
   return (gin);
}


double LnGamma (double x)
{
/* returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
   Stirling's formula is used for the central polynomial part of the procedure.


   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double f=0, fneg=0, z, lng;
   int nx=(int)x-1;


   if((double)nx==x && nx>0 && nx<10)
      lng=log((double)factorial(nx));
   else {
      if(x<=0) {
         error2("lnGamma not implemented for x<0");
         if((int)x-x==0) { puts("lnGamma undefined"); return(-1); }
         for (fneg=1; x<0; x++) fneg/=x;
         if(fneg<0) error2("strange!! check lngamma");
         fneg=log(fneg);
      }
      if (x<7) {
         f=1;  z=x-1;
         while (++z<7)  f*=z;
         x=z;   f=-log(f);
      }
      z = 1/(x*x);
      lng = fneg+ f + (x-0.5)*log(x) - x + .918938533204673 
             + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
                  +.083333333333333)/x;
   }
   return  lng;
}

/******************************************************************
                         tree functions
********************************************************************/
int ReadaTree (FILE *fTree,Tree *tree)
{
/* 
   Both names and numbers for species are accepted.  
   Species names are considered case-sensitive, with trailing blanks ignored.
*/
   	int cnode, cfather=-1, taxa = 0;  /* current node and father */
   	int inodeb=0;  /* node number that will have the next branch length */
   	int i, level=0, ch=' ';
   	char  skips[]="\"\'";
   	int nnode;   

	nnode = tree->ntaxa; 
   	FOR(i,2*(tree->ntaxa)-1) {
      		tree->nodes[i].father=-1;
      		tree->nodes[i].nson=0; 
		tree->nodes[i].sons[0] = -1;
		tree->nodes[i].sons[1] = -1;
   	}

   	while(ch != '(')
     	{
      		ch=fgetc(fTree);
     	}
   	ungetc(ch,fTree);

   	for (;;) {
      		ch = fgetc (fTree);
      		if (ch==EOF) return(-1);
      		else if (!isgraph(ch) || ch==skips[0] || ch==skips[1]) continue;
      		else if (ch=='(') {
         		level++;
         		cnode=nnode++;
   
         		if(nnode > 2*(tree->ntaxa)-1)
              		{
                  		printf("check tree: perhaps too many '('s");
                  		exit(-1);
              		}
         		if (cfather>=0) {
            			tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = cnode;
            			tree->nodes[cnode].father=cfather;
         		}
         		else
            			tree->root=cnode;
         		cfather=cnode;
      		}
      		else if (ch==')') { level--;  inodeb=cfather; cfather=tree->nodes[cfather].father; }
      		else if (ch==':') fscanf(fTree,"%lf",&tree->nodes[inodeb].brlens);
			else if (ch=='#') fscanf(fTree,"%lf",&tree->nodes[inodeb].lambda);
      		else if (ch==',') ;
      		else if (ch==';' && level!=0) 
         	{
            		printf("; in treefile");
            		exit(-1);
         	}
      		else if (isdigit(ch))
      		{ 
         		ungetc(ch, fTree); 
         		fscanf(fTree,"%d",&inodeb); 
         		inodeb--;
         		tree->nodes[inodeb].father=cfather;
         		tree->nodes[cfather].sons[tree->nodes[cfather].nson++]=inodeb;
      		}
		else if (isalpha(ch))
		{		
			i = 0;
			while(ch != ':')
			{
				if(ch != ' ')
					tree->nodes[taxa].taxaname[i++] = ch;
				ch = fgetc(fTree);
			}
			ungetc(ch, fTree);
			tree->nodes[taxa].father = cfather;
			tree->nodes[cfather].sons[tree->nodes[cfather].nson++] = taxa;
			inodeb = taxa;
			taxa++;
		}
      		if (level<=0) break;
   	}
   
   	for ( ; ; ) {
      		while(isspace(ch=fgetc(fTree)) && ch!=';' );
      		if (ch==':')       fscanf(fTree, "%lf", &tree->nodes[tree->root].brlens);
			else if (ch=='#')       fscanf(fTree, "%lf", &tree->nodes[tree->root].lambda);
      		else if (ch==';')  break;
      		else  { ungetc(ch,fTree);  break; }
   	}

   	if(nnode != 2*(tree->ntaxa)-1) 
	{ 
		printf(" # of nodes != %d\n",2*tree->ntaxa-1); 
		exit(-1);
	}
	return (0);
}


int PrintTree(Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted)
{

	char			*tempStr;
	int                     tempStrSize;

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
		{
		MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
		}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		{
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
		}

	SaveSprintf (&tempStr, &tempStrSize,"(");
	AddToPrintString (tempStr);
					
	WriteTreeToFile (tree, tree->root, showName, showBrlens, showTheta, isRooted);

	if(showTheta == YES) 
		SaveSprintf (&tempStr, &tempStrSize,")#%lf;\n",tree->nodes[tree->root].theta);
	else 
		SaveSprintf (&tempStr, &tempStrSize,");\n");
	AddToPrintString (tempStr);
	free (tempStr); 

	return (NO_ERROR);					

}

void WriteTreeToFile (Tree *tree, int inode, int showName, int showBrlens, int showTheta, int isRooted)
{

		char			*tempStr;
		int                      tempStrSize = 200;

		tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
		if (!tempStr)
			MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

		if (tree->nodes[inode].nson == 0)
			{
				if (showBrlens == YES)
				{
    				if(showName == 0) SaveSprintf (&tempStr, &tempStrSize, "%d:%lf", inode+1, tree->nodes[inode].brlens);
					else SaveSprintf (&tempStr, &tempStrSize, "%s:%lf", tree->nodes[inode].taxaname, tree->nodes[inode].brlens);
					AddToPrintString (tempStr);
					if((tree->nodes[inode].theta>0) && showTheta == YES) 
					{
						SaveSprintf (&tempStr, &tempStrSize, "#%lf", tree->nodes[inode].theta);
						AddToPrintString (tempStr);
					}
				}
				else
				{
					if(showName == 0) SaveSprintf (&tempStr, &tempStrSize, "%d", inode+1);
					else SaveSprintf (&tempStr, &tempStrSize, "%s", tree->nodes[inode].taxaname);
					AddToPrintString (tempStr);
				}
			}
		else
			{
				if (inode != tree->root)
				{
					SaveSprintf (&tempStr, &tempStrSize, "(");
					AddToPrintString (tempStr);
				}
				WriteTreeToFile (tree,tree->nodes[inode].sons[0], showName, showBrlens, showTheta, isRooted);
				SaveSprintf (&tempStr, &tempStrSize, ",");
				AddToPrintString (tempStr);
				WriteTreeToFile (tree,tree->nodes[inode].sons[1], showName, showBrlens, showTheta, isRooted);	
				if (inode != tree->root)
				{
					if (tree->nodes[inode].father == tree->root && isRooted == NO)
					{
						if (showBrlens == YES)
						{
							SaveSprintf (&tempStr, &tempStrSize, ",%d:%lf", tree->nodes[inode].father + 1, tree->nodes[tree->nodes[inode].father].brlens);
							AddToPrintString (tempStr);
			
							if((tree->nodes[tree->nodes[inode].father].theta>0) && showTheta == YES) 
							{
								SaveSprintf (&tempStr, &tempStrSize, "#%lf", tree->nodes[tree->nodes[inode].father].theta);
								AddToPrintString (tempStr);
							}
						}
						else
						{
							SaveSprintf (&tempStr, &tempStrSize, ",%d", tree->nodes[inode].father + 1);
							AddToPrintString (tempStr);
						}
					}
				
					if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
					{
						SaveSprintf (&tempStr, &tempStrSize,"):%lf", tree->nodes[inode].brlens);
						AddToPrintString (tempStr);
						if((tree->nodes[inode].theta > 0) && showTheta == YES)
						{
							SaveSprintf (&tempStr, &tempStrSize, "#%lf", tree->nodes[inode].theta);
							AddToPrintString (tempStr);
						}
					}
					else
					{
						SaveSprintf (&tempStr, &tempStrSize, ")");
						AddToPrintString (tempStr);
					}					
				}
			}
	free (tempStr);
		
		
}

int PrintNode(Tree *tree, int inode, int showBrlens, int showTheta, int isRooted)
{

	char			*tempStr;
	int                     tempStrSize;

	/* allocate the print string */
	printStringSize = 200;
	printString = (char *)SafeMalloc((size_t) (printStringSize * sizeof(char)));
	if (!printString)
		{
		MrBayesPrint ("%s   Problem allocating printString (%d)\n", spacer, printStringSize * sizeof(char));
		return (ERROR);
		}
	*printString = '\0';

	tempStrSize = 200;
	tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
	if (!tempStr)
		{
		MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));
		return (ERROR);
		}

	SaveSprintf (&tempStr, &tempStrSize,"(");
	AddToPrintString (tempStr);
					
	WriteNodeToFile (tree, tree->root, showBrlens, showTheta, isRooted);
	SaveSprintf (&tempStr, &tempStrSize,"):%d;",tree->root+1);
	AddToPrintString (tempStr);
	free (tempStr); 

	return (NO_ERROR);					

}

void WriteNodeToFile (Tree *tree, int inode, int showBrlens, int showTheta, int isRooted)
{

		char			*tempStr;
		int                      tempStrSize = 200;

		tempStr = (char *) SafeMalloc((size_t) (tempStrSize * sizeof(char)));
		if (!tempStr)
			MrBayesPrint ("%s   Problem allocating tempString (%d)\n", spacer, tempStrSize * sizeof(char));

		if (tree->nodes[inode].nson == 0)
			{
    				SaveSprintf (&tempStr, &tempStrSize, "%s:%d", tree->nodes[inode].taxaname,inode+1);
				AddToPrintString (tempStr);
			}
		else
			{
				if (inode != tree->root)
				{
					SaveSprintf (&tempStr, &tempStrSize, "(");
					AddToPrintString (tempStr);
				}
				WriteNodeToFile (tree,tree->nodes[inode].sons[0],  showBrlens, showTheta, isRooted);
				SaveSprintf (&tempStr, &tempStrSize, ",");
				AddToPrintString (tempStr);
				WriteNodeToFile (tree,tree->nodes[inode].sons[1], showBrlens, showTheta, isRooted);	
				if (inode != tree->root)
				{
					if (tree->nodes[inode].father == tree->root && isRooted == NO)
					{
						SaveSprintf (&tempStr, &tempStrSize, ",%d:%d", tree->nodes[inode].father + 1, tree->nodes[inode].father + 1);
						AddToPrintString (tempStr);
					}
				
					if (showBrlens == YES && isRooted == YES) /*tree->nodes[inode].father != tree->root)*/
					{
						SaveSprintf (&tempStr, &tempStrSize,"):%d", inode+1);
						AddToPrintString (tempStr);
					}
					else
					{
						SaveSprintf (&tempStr, &tempStrSize, ")");
						AddToPrintString (tempStr);
					}					
				}
			}
	free (tempStr);
		
		
}
void MrBayesPrint (char *format, ...)
{
	va_list                 ptr;

	va_start (ptr, format);

			vprintf (format, ptr);
			fflush (stdout);

	va_end(ptr);

}


int AddToPrintString (char *tempStr)

{

        size_t                  len1, len2;

        len1 = (int) strlen(printString);
        len2 = (int) strlen(tempStr);
        if (len1 + len2 + 5 > printStringSize)
                {
                printStringSize += len1 + len2 - printStringSize + 200;
                printString = realloc((void *)printString, printStringSize * sizeof(char));
                if (!printString)
                        {
                        MrBayesPrint ("%s   Problem reallocating printString (%d)\n", spacer, printStringSize * sizeof(char));
                        goto errorExit;
                        }
                }
        strcat(printString, tempStr);
#       if 0
        printf ("printString(%d) -> \"%s\"\n", printStringSize, printString);
#       endif
        return (NO_ERROR);

        errorExit:
                return (ERROR);

}

#define TARGETLENDELTA (100)

int SaveSprintf(char **target, int *targetLen, char *fmt, ...) {
  va_list    argp;
  int        len,retval;

  va_start(argp, fmt);
#ifdef VISUAL
  len = _vsnprintf(NULL, 0, fmt, argp);
#else
  len = vsnprintf(NULL, 0, fmt, argp);
#endif

  va_end(argp);

  if(len>*targetLen)
        {
/*        fprintf(stderr, "* increasing buffer from %d to %d bytes\n", *targetLen, len+TARGETLENDELTA); */
        *targetLen = len+TARGETLENDELTA; /* make it a little bigger */
            *target = realloc(*target, *targetLen);
        }

  va_start(argp,fmt);
  retval=vsprintf(*target, fmt, argp);
  va_end(argp);

/*   fprintf(stderr, "* savesprintf /%s/\n",*target); */
  return retval;
}

void *SafeMalloc(size_t s) 
{
        void *ptr = malloc(s);
        if(ptr==NULL)
                return NULL;
        return memset(ptr,0,s);
}


void PrintTreeToFile (FILE *file, Tree *tree)
{
	PrintTree(tree, tree->root, 1, 1, 1, 1);
	fprintf(file, "%s", printString);
	free (printString);
}

void PrintNodeToFile (FILE *file, Tree *tree)
{
	PrintNode(tree, tree->root, 1, 0, 1);
	fprintf(file, "%s", printString);
	free (printString);
}

void ShowNodes (Tree *tree)
{
	int i;

	for(i=0;i<(2*(tree->ntaxa)-1);i++)
		printf("node %d %d %d %d %f %d\n",i,tree->nodes[i].father,tree->nodes[i].sons[0],tree->nodes[i].sons[1], tree->nodes[i].brlens,tree->root);
	for(i=0;i<tree->ntaxa;i++) printf("%s ",tree->nodes[i].taxaname);
}

void findOffsprings (int *offsprings, Tree *tree, int inode)
{
	int son0, son1;

	if(inode < tree->ntaxa)
		offsprings[inode] = 1;
	else
	{
		son0 = tree->nodes[inode].sons[0];
		son1 = tree->nodes[inode].sons[1];
		findOffsprings (offsprings, tree, son0);
		findOffsprings (offsprings, tree, son1);
	}
}

double TreeHeight (Tree *tree)
{
	int father;
	double treeheight;

	father = tree->nodes[0].father;
	treeheight = tree->nodes[0].brlens;

	while(father != tree->root){
		treeheight += tree->nodes[father].brlens;
		father = tree->nodes[father].father;
	}
	return (treeheight);
}
