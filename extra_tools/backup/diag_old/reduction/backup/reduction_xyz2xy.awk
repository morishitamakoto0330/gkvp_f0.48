#!/bin/awk
# <Note>
#   Sum of 3D data over the third column.
#   The 3D data, wk[1:NR,1:NF], consists of x[1:nx], y[1:ny], z[1:nz], data[1:nx,1:ny,1:nz,1:nd] and comment or blank lines.
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
  if (NF>=4)                       nd = NF-3
  for ( ic=1; ic<=NF; ic++ ) {
    wk[NR,ic] = $ic
  }
}

END{
  if (ny!=-1) {
    nz = (NR-ncom)/((nx+1)*ny)
  } else {
    nz=1
    ny=(NR-ncom)/(nx+1)
  }

  ### Re-arrange the data ###
  for ( iz=1; iz<=nz; iz++ ) {
    for ( iy=1; iy<=ny; iy++ ) {
      for ( ix=1; ix<=nx; ix++ ) {
        for ( id=1; id<=nd; id++ ) {
          ir = ncom + (nx+1)*ny*(iz-1) + (nx+1)*(iy-1) + ix
          ic = id+3
          x[ix] = wk[ir,1]
          y[iy] = wk[ir,2]
          z[iz] = wk[ir,3]
          data[ix,iy,iz,id] = wk[ir,ic]
        }
      }
    }
  }

  ### Sum of the data over the third column ###
  for ( ix=1; ix<=nx; ix++ ) {
    for ( iy=1; iy<=ny; iy++ ) {
      for ( id=1; id<=nd; id++ ) {
        sum[ix,iy,id] = 0.0
      }
    }
  }
  for ( iz=1; iz<=nz; iz++ ) {
    for ( iy=1; iy<=ny; iy++ ) {
      for ( ix=1; ix<=nx; ix++ ) {
        for ( id=1; id<=nd; id++ ) {
          sum[ix,iy,id] = sum[ix,iy,id] + data[ix,iy,iz,id]
        }
      }
    }
  }

  printf( "%6s%d%6s%d%6s%d%6s%d\n", "# nx =", nx, ", ny =", ny, ", nz =", nz, ", nd =", nd )
  printf( "%17s%6s%17.7e%6s%17.7e\n", "#  Summed over z.", "From", z[1], "to", z[nz] )
  printf( "%17s%17s%17s\n", "#               x", "y", "data" )
  for ( iy=1; iy<=ny; iy++ ) {
    for ( ix=1; ix<=nx; ix++ ) {
      printf( "%17.7e%17.7e", x[ix], y[iy] )
      for ( id=1; id<=nd; id++ ) {
        printf( "%17.7e", sum[ix,iy,id] )
      #  printf( "%17.7e", sum[ix,iy,id]/nz ) # Average
      }
      printf( "\n" )
    }
    printf( "\n" )
  }
}
