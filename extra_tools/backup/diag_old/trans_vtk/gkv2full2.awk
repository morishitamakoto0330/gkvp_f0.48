#!/bin/awk
# <Note>
#   Transform of vtk file from GKV coordinate to Cartesian coordinate.
#
# <How to use>
#   awk -f gkv2cartesian.awk datafile > out.dat
#
#

BEGIN{
#  n_alp = 12
#  n1 = 41
#  n2 = 41
#  i3 = 0
}

(n1*n2*i3<=NR-12)&&(NR-12<n1*n2*(i3+1)){
  data[NR] = $0
}

END{
  for ( i_alp=0; i_alp<n_alp; i_alp++ ){
    for ( row=n1*n2*i3+12; row<n1*n2*(i3+1)+12; row++ ){
      print data[row]
    }
  }
}
