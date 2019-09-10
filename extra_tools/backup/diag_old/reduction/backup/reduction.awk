#!/bin/awk
# <Note>
#   Time derivative is calculated by 4th order central finite difference.
# <How to use>
#   awk -f reduction.awk filename > out.dat
#
# real :: time[0:ntime]
# real :: ky[0:global_ny]
# real :: data[0:ntime,0:global_ny,0:ndata]
#

BEGIN{
  global_ny=31
  time_sta=100.0
  time_end=180.0
  rows=0
}

($1!="#"&&$1!=""){
  it = int(rows/(global_ny+1))
  ntime = it
  time[it] = $1
  iy = rows%(global_ny+1)
  ky[iy] = $2
  ndata = NF-3
  for ( ic=0; ic<=ndata; ic++ ) {
    data[it,iy,ic] = $(ic+3)
  }
  rows++
}

END{
  for ( iy=0; iy<=global_ny; iy++ ) {
    for ( ic=0; ic<=ndata; ic++ ) {
      nave = 0
      average[iy,ic] = 0.0
    }
  }
  for ( it=0; it<=ntime; it++ ) {
    if ( time_sta<=time[it] && time[it]<=time_end ) {

      for ( iy=0; iy<=global_ny; iy++ ) {
        for ( ic=0; ic<=ndata; ic++ ) {
          average[iy,ic] = average[iy,ic] + data[it,iy,ic]
        }
      }
      nave++

    }
  }
  for ( iy=0; iy<=global_ny; iy++ ) {
    for ( ic=0; ic<=ndata; ic++ ) {
      average[iy,ic] = average[iy,ic]/nave
    }
  }

  printf( "%17s%17.7e%3s%17.7e\n",  \
          "# Averaged time =", time_sta, " ~ ", time_end )
  printf( "%17s%17s\n", "#              ky", "data" )
  for ( iy=0; iy<=global_ny; iy++ ) {
    printf( "%17.7e", ky[iy] )
    for ( ic=0; ic<=ndata; ic++ ) {
      printf( "%17.7e", average[iy,ic] )
    }
    printf( "\n" )
  }
}
