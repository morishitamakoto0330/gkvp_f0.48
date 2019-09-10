#!/bin/awk
# <Note>
#   Sum of 1D data.
#   The 1D data, wk[1:NR,1:NF], consists of x[1:nx], data[1:nx,1:nd] and comment or blank lines.
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment
#         x[ 1] data[ 1, 1] ... data[ 1,nd]         <- No blank line between comment and data.
#         x[ 2] data[ 2, 1] ... data[ 2,nd]
#           |  
#         x[nx] data[nx, 1] ... data[nx,nd]
#

BEGIN{
  ncom = 0   # Number of comment lines
  nx = -1    # Number of the first column
  nd = -1    # Number of data
}

{
  if ($1=="#") ncom++
  nx = NR-ncom
  if (NF>=2) nd = NF-1
  for ( ic=1; ic<=NF; ic++ ) {
    wk[NR,ic] = $ic
  }
}

END{
  ### Re-arrange the data ###
  for ( ix=1; ix<=nx; ix++ ) {
    for ( id=1; id<=nd; id++ ) {
      ir = ncom + ix
      ic = id+1
      x[ix] = wk[ir,1]
      data[ix,id] = wk[ir,ic]
    }
  }

  ### Sum of the data ###
  for ( id=1; id<=nd; id++ ) {
    sum[id] = 0.0
  }
  for ( ix=1; ix<=nx; ix++ ) {
    for ( id=1; id<=nd; id++ ) {
      sum[id] = sum[id] + data[ix,id]
    }
  }

  printf( "%6s%d%6s%d\n", "# nx =", nx, ", nd =", nd )
  printf( "%17s%6s%17.7e%6s%17.7e\n", "#  Summed over x.", "From", x[1], "to", x[nx] )
  printf( "%17s\n", "#            data" )
  for ( id=1; id<=nd; id++ ) {
    printf( "%17.7e", sum[id] )
  #  printf( "%17.7e", sum[id]/nx ) # Average
  }
  printf( "\n" )
}
