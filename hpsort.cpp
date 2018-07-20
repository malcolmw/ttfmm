#include <stdio.h>
#include <stdlib.h>
#include <deque>
#include "ttfmm.h"

void hpsort(unsigned long n, std::deque<GridNode> &ra)
{
  unsigned long i,ir,j,l;
  GridNode rra;

  if (n < 2) return;
  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1) {
      rra=ra[--l];
    } else {
      rra=ra[ir];
      ra[ir]=ra[1];
      if (--ir == 1) {
        ra[1]=rra;
        break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j].value < ra[j+1].value) j++;
      if (rra.value < ra[j].value) {
        ra[i]=ra[j];
        i=j;
        j <<= 1;
      } else break;
    }
    ra[i]=rra;
  }
}
