#!/bin/awk
# <Note>
#   Sort gkv.hhv.*** for gnuplot
# <How to use>
#   gawk -f sort_hhinvm.awk filename > out.dat
#

BEGIN{
  nprocv = 2
  nprocm = 2
  global_nv = 32
  global_nm = 31
  slngv = 2 * global_nv / nprocv
  slngm = ( global_nm + 1 ) / nprocm
}

NR<=2{
  coment[NR] = $0
}

NR>=3{
  tmp[NR-2] = $0
}

END{
  print coment[1]
  print coment[2]
  for( l = 0; l <= nprocm-1; l++ ){
    for( k = 0; k <= slngm-1; k++ ){
      for( j = 0; j <= nprocv-1; j++ ){
        for( i = 1; i <= slngv; i++ ){
          print tmp[i+slngv*slngm*j+slngv*k+slngv*slngm*nprocv*l]
        }
      }
      print ""
    }
  }
}
