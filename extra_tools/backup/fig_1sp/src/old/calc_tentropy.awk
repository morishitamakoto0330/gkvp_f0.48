#!/bin/awk
# <Note>
#   Time derivative is calculated by 4th order central finite difference.
# <How to use>
#   awk -f calc_entropybalance.awk filename > out.dat
#

BEGIN{
  dtout_eng = 0.1
}

{
   time[NR] = $1      # Time
      S[NR] = $2      # Entropy S
      W[NR] = $3      # Potential energy W
  ShQdt[NR] = $4      #
     hQ[NR] = $5      # Ion heat flux eta_i*Q
   SDdt[NR] = $6      #
      D[NR] = $7      # Collisional dissipation D
}

END{
  for(ir=3; ir<=NR-2; ir++){
    cef = 1 / ( 12 * dtout_eng )
    dSdt[ir] = cef * ( - S[ir+2] + 8 * S[ir+1] - 8 * S[ir-1] + S[ir-2] )
    dWdt[ir] = cef * ( - W[ir+2] + 8 * W[ir+1] - 8 * W[ir-1] + W[ir-2] )
  }

  printf( "%17s%17s%17s%17s%17s%17s%17s%17s%17s\n",  \
          "#            time", "S", "W", "ShQdt", "SDdt", "dSdt", "dWdt", "hQ", "D" )
  for(ir=1; ir<=2; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17s%17s%17.7e%17.7e\n",  \
            time[ir], S[ir], W[ir], ShQdt[ir], SDdt[ir],             \
                      "NaN", "Nan", hQ[ir], D[ir]                   )
  }
  for(ir=3; ir<=NR-2; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            time[ir], S[ir], W[ir], ShQdt[ir], SDdt[ir],                 \
                      dSdt[ir], dWdt[ir], hQ[ir], D[ir]                 )
  }
  for(ir=NR-1; ir<=NR; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17s%17s%17.7e%17.7e\n",  \
            time[ir], S[ir], W[ir], ShQdt[ir], SDdt[ir],             \
                      "NaN", "Nan", hQ[ir], D[ir]                   )
  }
}
