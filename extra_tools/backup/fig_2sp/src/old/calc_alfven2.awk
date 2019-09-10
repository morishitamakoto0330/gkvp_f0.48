#!/bin/awk
# <Note>
#   Calculate the frequency and growthrate of a standing wave
# <How to use>
#   awk -f calc_alfvenwave.awk filename > out.dat
#

BEGIN{
  t_rec = 11.5        # recurrence time
  dcolumn = 0        # number of dumped columns
  tcolumn = 4        # number of trimmed columns
  drow = 1           # number of dumped rows
  trow = 1000        # number of trimmed columns
}

(drow+1) <= NR && NR <= (drow+trow) {
  data[NR-drow,1] = $(dcolumn+1)       # ic=1 is time.
  for ( ic = 2; ic <= tcolumn; ic++ ) {
    data[NR-drow,ic] = $(dcolumn+ic)   # ic=2,3,4 are S, WE, WM
  }
}

END{
  pi = 3.1415926535897
  prev = data[1,3]
  trig = 0
  flag = 0
  for ( ir = 1; ir <= trow; ir++ ) {
    if ( data[ir,3] > prev ) {
      if ( trig == -1 ) {
        flag++
        t_0[flag] = data[ir,1]
      }
      trig = +1
    } else {
      trig = -1
    }
    prev = data[ir,3]
    norm[flag] = norm[flag] + data[ir,3]
    if ( data[ir,1] > t_rec ) break
  }
#  for ( ir = 1; ir <= flag; ir++ ) {
#    printf( "%d%17.7e%17.7e\n", ir, t_0[ir], norm[ir] )
#  }
  period = 2 * ( t_0[flag] - t_0[1] ) / ( flag - 1 )
  freq = 2 * pi / period
  growth = ( log(norm[flag-1]) - log(norm[1]) ) / ( ( flag - 2 ) * period )
  printf( "%17.7e%17.7e\n", freq, growth )
}
