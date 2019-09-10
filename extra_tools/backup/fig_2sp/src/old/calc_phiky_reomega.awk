#!/bin/awk
# <Note>
#   Derivative is calculated by least square method.
# <How to use>
#   awk -f growthrate.awk filename > out.dat
#

BEGIN{
  kymin = 0.07
  dcolumn=0         # number of dumped columns
  tcolumn=9        # number of trimmed columns
  drow=0            # number of dumped rows
  trow=10           # number of trimmed columns
}

(drow+1)<=NR && NR<=(drow+trow){
  data[NR-drow,1]=$(dcolumn+1)          # ic=1 is time.
  for(ic=2; ic<=tcolumn; ic++){
    data[NR-drow,ic]=$(dcolumn+ic)      # phase=-reomega*t+phase0
  }
}

END{
  for(ic=2; ic<=tcolumn; ic++){
    B=0
    C=0
    D=0
    E=0
    for(ir=1; ir<=trow; ir++){
      B=B+data[ir,1]**2             # sum_i(xi**2)
      C=C+data[ir,ic]               # sum_i(yi)
      D=D+data[ir,1]*data[ir,ic]    # sum_i(xi*yi)
      E=E+data[ir,1]                # sum_i(xi)
    }
    a=(trow*D-C*E)/(trow*B-E**2)    # function form is 'y=a*x+y0'
    y0=(B*C-D*E)/(trow*B-E**2)

    ky=kymin*(ic-1)                  # ky = 2pi/Ly * iky
    reomega=-a                      # reomega = -a
    printf("%17.7e%17.7e\n",ky,reomega)
  }
}
