#!/bin/awk
# <Note>
#   Transform of vtk file from GKV coordinate to Cartesian coordinate.
#
# <How to use>
#   awk -f gkv2cartesian.awk datafile > out.dat
#
#

BEGIN{
#  r_minor = 84.0       # Minor radius [in rho unit]
#  eps_r   = 0.18       # Inverse aspect ratio (=r/R_0)
#  q_0     = 1.4        # Safety factor
#  s_hat   = 0.78       # Magnetic shear
#  n_alp   = 12
}

(NR==1){ print $0 }
(NR==2){
  print $0
  lx = $2
  ly = $4
  lz = $6
}
(NR==3){ print $0 }
(NR==4){ printf "%s\n", "DATASET STRUCTURED_GRID" }
(NR==5){
  print $1, 1, $3*n_alp, $4
  nxw       = ($2-1)/2
  nyw       = ($3-1)/2
  global_nz = ($4-1)/2
  dx = lx / nxw
  dy = ly / nyw
  dz = lz / global_nz
  q_bar = q_0 * sqrt( 1.0 - eps_r*eps_r )
  pi = 3.141592653589793

  printf("%s%17d%s\n", "POINTS ", 1*(2*nyw+1)*n_alp*(2*global_nz+1), " float" )
  for ( iz=-global_nz; iz<=global_nz; iz++ ){
    zz = dz * iz
    theta = 2.0 * atan( sqrt( (1.0+eps_r)/(1.0-eps_r) ) * tan( 0.5*zz ) )
    for ( i_alp=0; i_alp<n_alp; i_alp++ ){
      for ( my=0; my<=2*nyw; my++ ){
        yy = -ly + dy * my
        #for ( mx=0; mx<=2*nxw; mx++ ){
	mx = 0
          xx = -lx + dx * mx
          sr = r_minor + xx
          mr = r_minor / eps_r + sr * cos( theta )
          zeta = q_bar * ( zz + ( s_hat * xx * zz - yy )/r_minor ) \
                                        - 2.0 * pi * i_alp / n_alp
          x_car = mr * cos( zeta )
          y_car = mr * sin( zeta )
          z_car = sr * sin( theta )
          printf("%17.7e%17.7e%17.7e\n",x_car,y_car,z_car)
        #}
      }
    }
  }
}
(NR==9){ print $1, $2*n_alp/(2*nxw+1) }
(NR==10){ print $0 }
(NR==11){ print $0 }

function tan(x) {
  return sin(x)/cos(x)
}
function atan(x) {
  return atan2(x,1)
}
