/*                                                                              */
/*      WeICME (Wedge Integrated Computational Materials Engineering)           */
/*                 - A 3-dimensional finite element program.                    */   
/*     Developed and maintained by Shenzhen Wedge Central                       */
/*    South Research Institute co., Ltd., Shenzhen, China                       */   
/*     Copy Right 2019-2023.                                                      */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "WeICME.h"

void frdgeneralvector(double *v,ITG *iset,ITG *ntrans,char * filabl,
               ITG *nkcoords,
               ITG *inum,char *m1,ITG *inotr,double *trab,double *co,
               ITG *istartset,ITG *iendset,ITG *ialset,ITG *mi,ITG *ngraph,
               FILE *f1,char *output,char *m3){

    ITG i,j,k,l,m,nksegment;
      
  int iw;

  float ifl;
  
  double a[9];

  if(*iset==0){
    if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)){
      for(i=0;i<*nkcoords;i++){
	if(inum[i]<=0) continue;
	if(strcmp1(output,"asc")==0){
	    if(mi[1]==4){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
	    }else if(mi[1]==5){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5]);
	    }else if(mi[1]==6){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
	    }
	}else{
	  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	  for(j=1;j<=mi[1];j++){
	      ifl=(float)v[(mi[1]+1)*i+j];fwrite(&ifl,sizeof(float),1,f1);
	  }
	}
      }
    }else{
      for(i=0;i<*nkcoords;i++){
	if(inum[i]<=0) continue;
	if(inotr[2*i]==0){
	  if(strcmp1(output,"asc")==0){
	    if(mi[1]==4){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
	    }else if(mi[1]==5){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5]);
	    }else if(mi[1]==6){
		fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			(float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
	    }
	  }else{
	    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
	    for(j=1;j<=mi[1];j++){
		ifl=(float)v[(mi[1]+1)*i+j];fwrite(&ifl,sizeof(float),1,f1);
	    }
	  }
	}else{
	    printf("*WARNING in frdgeneralvector:\n");
	    printf("         no output in local coordinates allowed\n");
	    printf("         for generalized vectors\n");
	    printf("         output request ist not performed;\n");
	}
      }
    }
  }else{
    nksegment=(*nkcoords)/(*ngraph);
    for(k=istartset[*iset-1]-1;k<iendset[*iset-1];k++){
      if(ialset[k]>0){
	for(l=0;l<*ngraph;l++){
	  i=ialset[k]+l*nksegment-1;
	  if(inum[i]<=0) continue;
	  if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)||(inotr[2*i]==0)){
	      if(strcmp1(output,"asc")==0){
		  if(mi[1]==4){
		      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
		  }else if(mi[1]==5){
		      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			      (float)v[(mi[1]+1)*i+5]);
		  }else if(mi[1]==6){
		      fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
			      (float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
			      (float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
			      (float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
		  }
	      }else{
		  iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
		  for(j=1;j<=mi[1];j++){
		      ifl=(float)v[(mi[1]+1)*i+j];fwrite(&ifl,sizeof(float),1,f1);
		  }
	      }
	  }else{
	      printf("*WARNING in frdgeneralvector:\n");
	      printf("         no output in local coordinates allowed\n");
	      printf("         for generalized vectors\n");
	      printf("         output request ist not performed;\n");
	  }
	}
      }else{
	l=ialset[k-2];
	do{
	  l-=ialset[k];
	  if(l>=ialset[k-1]) break;
	  for(m=0;m<*ngraph;m++){
	    i=l+m*nksegment-1;
	    if(inum[i]<=0) continue;
	    if((*ntrans==0)||(strcmp1(&filabl[5],"G")==0)||(inotr[2*i]==0)){
		if(strcmp1(output,"asc")==0){
		    if(mi[1]==4){
			fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
				(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
				(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4]);
		    }else if(mi[1]==5){
			fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
				(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
				(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
				(float)v[(mi[1]+1)*i+5]);
		    }else if(mi[1]==6){
			fprintf(f1,"%3s%10" ITGFORMAT "%12.5E%12.5E%12.5E%12.5E%12.5E%12.5E\n",m1,i+1,
				(float)v[(mi[1]+1)*i+1],(float)v[(mi[1]+1)*i+2],
				(float)v[(mi[1]+1)*i+3],(float)v[(mi[1]+1)*i+4],
				(float)v[(mi[1]+1)*i+5],(float)v[(mi[1]+1)*i+6]);
		    }
		}else{
		    iw=(int)(i+1);fwrite(&iw,sizeof(int),1,f1);
		    for(j=1;j<=mi[1];j++){
			ifl=(float)v[(mi[1]+1)*i+j];fwrite(&ifl,sizeof(float),1,f1);
		    }
		}
	    }else{
	      printf("*WARNING in frdgeneralvector:\n");
	      printf("         no output in local coordinates allowed\n");
	      printf("         for generalized vectors\n");
	      printf("         output request ist not performed;\n");
	    }
	  }
	}while(1);
      }
    }
  }
      
  if(strcmp1(output,"asc")==0)fprintf(f1,"%3s\n",m3);

  return;

}

