#!/bin/awk
# <Note>
#   Create gnuplot file
#
# <Data format>
#   z,time,data;  ky,time,data
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
#  time_spacing = 1
  ncom = -1         # Number of comment lines
  nx = -1           # Number of the first column
  ny = -1           # Number of the second column
  nd = -1           # Number of data
}

{
  if (ncom==-1&&$1!="#"&&$1!="") { ncom=NR-1; wx=$1; wy=$2 }
  if (nx==-1&&$1==wx)              nx = NR-ncom-2
#  if (NF>=3)                       nd = NF-2
#  for ( ic=1; ic<=NF; ic++ ) {
#    wk[NR,ic] = $ic
    wk[NR,2] = $2 #Revised
#  }
}

END{
  if (nx!=-1) {
    ny = (NR-ncom)/(nx+1)
  } else {
    ny = 1
    nx = NR-ncom-1
  }

  ### Re-arrange the data ###
  for ( iy=1; iy<=ny; iy=iy+time_spacing ) {
#    for ( ix=1; ix<=nx; ix++ ) {
     ix=1 #Revised
#      for ( id=1; id<=nd; id++ ) {
        ir = ncom + (nx+1)*(iy-1) + ix
#        ic = id+2
#        x[ix] = wk[ir,1]
        y[iy] = wk[ir,2]
#        data[ix,iy,id] = wk[ir,ic]
#      }
#    }
  }

  for ( iy=1; iy<=ny; iy=iy+time_spacing ) {
    ll = (nx+1)
    le = ncom + (nx+1)*iy
    printf( "%s%s%s%08d%s\n", "# set output '", FILENAME, "_l", iy-1, ".jpg'" )
    printf( "%s%17.7f%s\n", "  set title 'time =", y[iy], "'" )
    printf( "%s%d%s%s%s%d%s\n", "  plot '< head -", le, " ", FILENAME, " | tail -", ll, "' u 1:3 w lp, '' u 1:4 w lp" )
    printf( "%s\n", "  pause 0.2" )
  }
}
