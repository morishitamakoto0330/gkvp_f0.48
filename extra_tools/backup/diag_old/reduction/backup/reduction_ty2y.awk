#!/bin/awk
# <Note>
#   Average 2D data over the first column.
#   The 2D data, wk[1:NR,1:NF], consists of x1[1:n1], x2[1:n2], data[1:n1,1:n2,1:nd] and commented-out or blank lines.
# <How to use>
#   sed -n "ir_sta,ir_endp" filename | awk -f reduction.awk > out.dat
#

BEGIN{
  ncom = 0   # Number of commented-out lines
  n1 = -1    # Number of the second column
  n2 = -1    # Number of the first column
  nd = -1    # Number of data
}

{
  if ($1=="#") ncom++
  if ($1==""&&n2==-1) n2 = NR-ncom-1
  for ( ic=1; ic<=NF; ic++ ) {
    wk[NR,ic] = $ic
  }
  if (n2!=-1) n1=(NR-ncom)/(n2+1)
  if (NF>=3) nd = NF-2
}

END{
  ### Re-arrange the data ###
  for ( i1=1; i1<=n1; i1++ ) {
    for ( i2=1; i2<=n2; i2++ ) {
      for ( id=1; id<=nd; id++ ) {
        ir = ncom + (i1-1)*(n2+1) + i2 # where "1" added to n2 corresponds to blank line
        ic = id+2
        x1[i1] = wk[ir,1]
        x2[i2] = wk[ir,2]
        data[i1,i2,id] = wk[ir,ic]
      }
    }
  }

  ### Average the data over the first column ###
  nave = 0
  for ( i2=1; i2<=n2; i2++ ) {
    for ( id=1; id<=nd; id++ ) {
      ave[i2,id] = 0.0
    }
  }
  for ( i1=1; i1<=n1; i1++ ) {
    for ( i2=1; i2<=n2; i2++ ) {
      for ( id=1; id<=nd; id++ ) {
        ave[i2,id] = ave[i2,id] + data[i1,i2,id]
      }
    }
    nave++
  }
  for ( i2=1; i2<=n2; i2++ ) {
    for ( id=1; id<=nd; id++ ) {
      ave[i2,id] = ave[i2,id]/nave
    }
  }

  printf( "%17s%17.7e%4s%17.7e\n", "#  Averaged x1 = ", x1[1], " to ", x1[n1] )
  printf( "%17s%17s\n", "#              x2", "data" )
  for ( i2=1; i2<=n2; i2++ ) {
    printf( "%17.7e", x2[i2] )
    for ( id=1; id<=nd; id++ ) {
      printf( "%17.7e", ave[i2,id] )
    }
    printf( "\n" )
  }
}
