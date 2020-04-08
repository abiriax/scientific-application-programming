/*
In the name of Allah, the Lord of Mercy in this life and the next
program title: TheoreticalSolidStateBscProject.C
author: Mohammed Basith Awan
date: 21/11/2010
objective: program is designed to investigate the nature of the patterns of
hexagonal (honeycomb) crystalline structures using Baxter's 3 colour theorem
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "strings.h"

#define L 3  /* linear system size */

/* random number initialising and call functions (see ran_mar.c) */
void rmarin(int ij,int kl,double *ran_pard,int *ran_pari);
void ranmar(double *rvec,int len,double *ran_pard,int *ran_pari);

/* checks that the 3-colour constraint is respected across the system */
void MatrixCheck(int mat[L][L][3]);

int main(int argc, char **argv)
{
FILE *fid;
int ind_hist;
int seed1,seed2;
int ii,jj,aa;
int ii0,jj0,col1,col2,aa1,aa2;
double ran_pard_tmp,ran_pard[101];
int ran_pari_tmp,ran_pari[3];
double rvec[4];
int mat[L][L][3];
int mat_store[L][L][3];
int ind;


/********************* Random number initialisation ***********************/
ind_hist=1;
fid=fopen("DATA/rand_seed50k.dat","r");
for (ii=0;ii<ind_hist;ii++)
  fscanf(fid,"%i %i\n",&seed1,&seed2);
fclose(fid);
rmarin(seed1,seed2,ran_pard,ran_pari);

printf("first step!\n");


/********************** Initialization of the system ***********************/
for (ii=0;ii<L;ii++)
  for (jj=0;jj<L;jj++)
    for (aa=0;aa<3;aa++)
      mat[ii][jj][aa]=aa;

printf("second step!\n");


/************************* 3 colour contraint check ************************/
MatrixCheck(mat);
printf("third step!\n");


/***************************************************************************/
/********************** BEGINNING of Loop Iteration ************************/
/***************************************************************************/
for (ind=0;ind<100;ind++)
{
printf("ind = %i\n",ind);

/***************************** Loop Update *********************************/
for (ii=0;ii<L;ii++)
  for (jj=0;jj<L;jj++)
    for (aa=0;aa<3;aa++)
      mat_store[ii][jj][aa]=mat[ii][jj][aa];

ranmar(rvec,4,ran_pard,ran_pari);
ii0=(int) (rvec[0]*L);
jj0=(int) (rvec[1]*L);
col1=(int) (rvec[2]*3.0);
col2=(col1+((int) (rvec[3]*2.0))+1)%3;

for (aa=0;aa<3;aa++)
  {
  if (mat[ii0][jj0][aa]==col1)
    aa1=aa;
  else (mat[ii0][jj0][aa]==col2)
    aa2=aa;
  }

ii=ii0;
jj=jj0;
aa=aa1;

do
  {
  mat[ii][jj][aa]=col2;
  switch (aa)
    {
    case 0: /* ii,jj, i.e. vertex 0 is col1 */
    if (mat[ii][(jj-1+L)%L][1]==col2) /* ii,jj-1,1 is col2 */
        {     
        ii=ii;
        jj=(jj-1+L)%L;
        mat[ii][jj][1]=col1;
        	if (mat[ii][jj][0]==col1)
        	{
        	 aa=0;	
        	}
        	else
        	{
      		 aa=2;	
        	}
        }
        
    else /* ii+1,jj-1,2 is col2 */
   		  {
          ii=(ii+1)%L;
          jj=(jj-1+L)%L;
          mat[ii][jj][2]=col1;
          if (mat[ii][jj][0]==col1)
          {
           aa=0;	
          }
          
          else
          {
  	       aa=1;
          }
		      
    }
      break;
      
    case 1: //ii,jj-1, i.e. vertex 1 is col2
	if (mat[ii][(jj+L+5)%L][1]==col) /* {fill in} is col2 */
		{
		ii=ii;
		jj	
		}
	else
		{
		if
		{
		 aa = 0;
		}
		
		else
		{
		 aa = 2;
		}
		}
      break;
      
 	if (mat[ii][(jj)%][1]==col) /* ii,jj, is col2 */
 	
    case 2://ii+1,jj-1, i.e. vertex 2  is col2
    if (mat[ii][(jj+L+5)%L][1]==col) /* ii,jj+1, is col2 */
      break;
    default:
      printf("There is a problem!\n");
      break;
    }
  }
while ( (ii != ii0) || (jj != jj0) );


/************************* 3 colour contraint check ************************/
MatrixCheck(mat);
MatrixCheck(mat_store);

}
/***************************************************************************/
/************************* END of Loop Iteration ***************************/
/***************************************************************************/

/***************************** Matrix storage ******************************/
fid=fopen("DATA/mat.dat","w");
for (ii=0;ii<L;ii++)
  for (jj=0;jj<L;jj++)
    fprintf(fid,"%i %i %i\n",mat[ii][jj][0],mat[ii][jj][1],mat[ii][jj][2]);
fclose(fid);



//printf("%lf %lf %lf\n",rvec[0],rvec[1],rvec[2]);

}
