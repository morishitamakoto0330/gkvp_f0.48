#!/bin/awk
# <Note>
#   Sum of 2D data over the first column.
#   The 2D data, wk[1:NR,1:NF], consists of x[1:nx], y[1:ny], data[1:nx,1:ny,1:nd] and comment or blank lines.
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment
#         x[ 1] y[ 1] data[ 1, 1, 1] ... data[ 1, 1,nd]          <- No blank line between comment and data.
#         x[ 2] y[ 1] data[ 2, 1, 1] ... data[ 2, 1,nd]
#           |     |
#         x[nx] y[ 1] data[nx, 1, 1] ... data[nx, 1,nd]
#                                                                <- A blank line between data blocks.
#         x[ 1] y[ 2] data[ 1, 2, 1] ... data[ 1, 2,nd]
#         x[ 2] y[ 2] data[ 2, 2, 1] ... data[ 2, 2,nd]
#           |     |
#         x[nx] y[ 2] data[nx, 2, 1] ... data[nx, 2,nd]
#
#                          |
#
#         x[ 1] y[ny] data[ 1,ny, 1] ... data[ 1,ny,nd]
#         x[ 2] y[ny] data[ 2,ny, 1] ... data[ 2,ny,nd]
#           |     |
#         x[nx] y[ny] data[nx,ny, 1] ... data[nx,ny,nd]
#                                                                <- A blank line at the end.
#

BEGIN{
  ncom = -1         # Number of comment lines
  nx = -1           # Number of the first column
  ny = -1           # Number of the second column
  nd = -1           # Number of data
}

{
  if (ncom==-1&&$1!="#"&&$1!="") { ncom=NR-1; wx=$1; wy=$2 }
  if (nx==-1&&$1==wx)              nx = NR-ncom-2
  if (NF>=3)                       nd = NF-2
  for ( ic=1; ic<=NF; ic++ ) {
    wk[NR,ic] = $ic
  }
}

END{
  if (nx!=-1) {
    ny = (NR-ncom)/(nx+1)
  } else {
    ny = 1
    nx = NR-ncom-1
  }

  ### Re-arrange the data ###
  for ( iy=1; iy<=ny; iy++ ) {
    for ( ix=1; ix<=nx; ix++ ) {
      for ( id=1; id<=nd; id++ ) {
        ir = ncom + (nx+1)*(iy-1) + ix
        ic = id+2
        x[ix] = wk[ir,1]
        y[iy] = wk[ir,2]
        data[ix,iy,id] = wk[ir,ic]
      }
    }
  }

  ### Sum of the data over the first column ###
  for ( iy=1; iy<=ny; iy++ ) {
    for ( id=1; id<=nd; id++ ) {
      sum[iy,id] = 0.0
    }
  }
  for ( iy=1; iy<=ny; iy++ ) {
    for ( ix=1; ix<=nx; ix++ ) {
      for ( id=1; id<=nd; id++ ) {
        sum[iy,id] = sum[iy,id] + data[ix,iy,id]
      }
    }
  }

  printf( "%6s%d%6s%d%6s%d\n", "# nx =", nx, ", ny =", ny, ", nd =", nd )
  printf( "%17s%6s%17.7e%6s%17.7e\n", "#  Summed over x.", "From", x[1], "to", x[nx] )
  printf( "%17s%17s\n", "#               y", "data" )
  for ( iy=1; iy<=ny; iy++ ) {
    printf( "%17.7e", y[iy] )
    for ( id=1; id<=nd; id++ ) {
      printf( "%17.7e", sum[iy,id] )
    #  printf( "%17.7e", sum[iy,id]/nx ) # Average
    }
    printf( "\n" )
  }
}
