#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef int IDL_LONG; /* #include "IDL_export.h" */
void hunt(float *xx, IDL_LONG n, float x, int *jlo)
{
  int jm,jhi,inc,ascnd,n1;

  n1=n-1;
  ascnd=(xx[n1] >= xx[0]);
  if (*jlo < 0 || *jlo >= n) {     /* no useful input -> bisection */
    *jlo=0;
    jhi= n;
  } else {
    inc=1;                           /* hunting increment */
    if ((x >= xx[*jlo]) == ascnd) {     /* hunt up */
      if (*jlo == n1) return;
      jhi = (*jlo)+1;
      while ((x > xx[jhi]) == ascnd) { /* not done yet */
	inc+=inc;                    /* double increment */
	jhi=(*jlo)+inc;
	if (jhi >n1) { 
	  jhi = n;                 /* Done hunting since off end of table */
	  break;
	}                            /* try again */
      }                    /* done hunting, value bracketed */
    } else {
      if (*jlo == 0) {               /* hunt down */
	return;
      }
      jhi = (*jlo)--;
      while((x < xx[*jlo]) == ascnd){     /* not done yet */
	jhi=(*jlo);
	inc*=2;
	if (inc >= jhi) {          /* off end of table */
	  *jlo=0;
	  break;
	}
	else *jlo = jhi-inc;
      }                           /* try again */
    }                             /* hunt done */
  }                               /* hunt done, begin final bisection */
  while( (jhi-(*jlo)) != 1) {
    jm=(jhi+(*jlo))/2;
    if ((x > xx[jm]) == ascnd)
      *jlo = jm;
    else
      jhi = jm;
  }
}

