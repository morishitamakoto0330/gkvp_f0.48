#!/bin/awk
# <Note>
#   Create gnuplot file
#
# <Data format>
#   x,y,time,data;  kx,ky,time,data
#
# <How to use>
#   awk -f reduction.awk datafile > out.dat
#
# <Example of datafile>
#       # Comment
#         x[ 1] y[ 1] z[ 1] data[ 1, 1, 1, 1] ... data[ 1, 1, 1,nd]         <- No blank line between comment and data.
#           |     |     |   
#         x[nx] y[ 1] z[ 1] data[nx, 1, 1, 1] ... data[nx, 1, 1,nd]
#                                                                           <- A blank line between data blocks.
#                                |
#         x[ 1] y[ny] z[ 1] data[ 1,ny, 1, 1] ... data[ 1,ny, 1,nd]
#           |     |     |   
#         x[nx] y[ny] z[ 1] data[nx,ny, 1, 1] ... data[nx,ny, 1,nd]
#
#                                |
#         x[ 1] y[ 1] z[nz] data[ 1, 1,nz, 1] ... data[ 1, 1,nz,nd]
#           |     |     |   
#         x[nx] y[ 1] z[nz] data[nx, 1,nz, 1] ... data[nx, 1,nz,nd]
#
#                                |
#         x[ 1] y[ny] z[nz] data[ 1,ny,nz, 1] ... data[ 1,ny,nz,nd]
#           |     |     |   
#         x[nx] y[ny] z[nz] data[nx,ny,nz, 1] ... data[nx,ny,nz,nd]
#                                                                           <- A blank line at the end.
#


BEGIN{
#  time_spacing = 1
  ncom = -1         # Number of comment lines
  nx = -1           # Number of the first column
  ny = -1           # Number of the second column
  nz = -1           # Number of the third column
  nd = -1           # Number of data
}

{
  if (ncom==-1&&$1!="#"&&$1!="") { ncom=NR-1; wx=$1; wy=$2; wz=$3 }
  if (nx==-1&&$1==wx)              nx = NR-ncom-2
  if (ny==-1&&nx!=-1&&$2==wy)      ny = (NR-ncom-1)/(nx+1)
#  if (NF>=4)                       nd = NF-3
#  for ( ic=1; ic<=NF; ic++ ) {
#    wk[NR,ic] = $ic
    wk[NR,3] = $3 #Revised
#  }
}

END{
  if (ny!=-1) {
    nz = (NR-ncom)/((nx+1)*ny)
  } else {
    nz=1
    ny=(NR-ncom)/(nx+1)
  }

  ### Re-arrange the data ###
  for ( iz=1; iz<=nz; iz=iz+time_spacing ) {
#    for ( iy=1; iy<=ny; iy++ ) {
#      for ( ix=1; ix<=nx; ix++ ) {
     iy=1; ix=1 #Revised
#        for ( id=1; id<=nd; id++ ) {
          ir = ncom + (nx+1)*ny*(iz-1) + (nx+1)*(iy-1) + ix
#          ic = id+3
#          x[ix] = wk[ir,1]
#          y[iy] = wk[ir,2]
          z[iz] = wk[ir,3]
#          data[ix,iy,iz,id] = wk[ir,ic]
#        }
#      }
#    }
  }

  for ( iz=1; iz<=nz; iz=iz+time_spacing ) {
    ll = (nx+1)*ny
    le = ncom + (nx+1)*ny*iz
    printf( "%s%s%s%08d%s\n", "# set output '", FILENAME, "_l", iz-1, ".jpg'" )
    printf( "%s%17.7f%s\n", "  set title 'time =", z[iz], "'" )
    printf( "%s%d%s%s%s%d%s\n", "  splot '< head -", le, " ", FILENAME, " | tail -", ll, "' u 1:2:4" )
    printf( "%s\n", "  pause 0.2" )
  }
}
