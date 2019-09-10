#!/bin/awk
# <Note>
#   Time derivative is calculated by 4th order central finite difference.
# <How to use>
#   awk -f calc_entropybalance.awk filename > out.dat
#

BEGIN{
  R0_Lt =
}

{
      t[NR] = $1                              # Time
  Ss_nz[NR] = $2                              #      Gyrocenter entropy S_s
  Ss_zf[NR] = $3                              #
  WE_nz[NR] = $4                              #         Electric energy W_E
  WE_zf[NR] = $5                              #
  RE_nz[NR] = $6                              #    Electric interaction R_sE
  RE_zf[NR] = $7                              #
  NE_nz[NR] = $8                              #   Nonlinear interaction N_sE
  NE_zf[NR] = $9                              #
  Ds_nz[NR] = $10                             # Collisional dissipation D_s
  Ds_zf[NR] = $11                             #
     QE[NR] = R0_Lt * $12                     #      Electric heat flux Q_sE
}

END{
  for(ir=3; ir<=NR-2; ir++){
### for uniform time steps ###
#    dtout_eng =
#    cef = 1 / ( 12 * dtout_eng )
#    dSsdt[ir] = cef * ( - Ss[ir+2] + 8 * Ss[ir+1] - 8 * Ss[ir-1] + Ss[ir-2] )
#    dWEdt[ir] = cef * ( - WE[ir+2] + 8 * WE[ir+1] - 8 * WE[ir-1] + WE[ir-2] )
### for non-uniform time steps ###
    cefm2 = ( (t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+1])*(t[ir  ]-t[ir+2]) )                    \
          / ( (t[ir-2]-t[ir-1])*(t[ir-2]-t[ir  ])*(t[ir-2]-t[ir+1])*(t[ir-2]-t[ir+2]) )
    cefm1 = ( (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir+1])*(t[ir  ]-t[ir+2]) )                    \
          / ( (t[ir-1]-t[ir-2])*(t[ir-1]-t[ir  ])*(t[ir-1]-t[ir+1])*(t[ir-1]-t[ir+2]) )
    cefp0 = (   (t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+1])*(t[ir  ]-t[ir+2])                    \
              + (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir+1])*(t[ir  ]-t[ir+2])                    \
              + (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+2])                    \
              + (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+1]) )                  \
          / ( (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+1])*(t[ir  ]-t[ir+2]) )
    cefp1 = ( (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+2]) )                    \
          / ( (t[ir+1]-t[ir-2])*(t[ir+1]-t[ir-1])*(t[ir+1]-t[ir  ])*(t[ir+1]-t[ir+2]) )
    cefp2 = ( (t[ir  ]-t[ir-2])*(t[ir  ]-t[ir-1])*(t[ir  ]-t[ir+1]) )                    \
          / ( (t[ir+2]-t[ir-2])*(t[ir+2]-t[ir-1])*(t[ir+2]-t[ir  ])*(t[ir+2]-t[ir+1]) )

    dSsdt_nz[ir] = cefm2*Ss_nz[ir-2] + cefm1*Ss_nz[ir-1] + cefp0*Ss_nz[ir  ] + cefp1*Ss_nz[ir+1] + cefp2*Ss_nz[ir+2]
    dSsdt_zf[ir] = cefm2*Ss_zf[ir-2] + cefm1*Ss_zf[ir-1] + cefp0*Ss_zf[ir  ] + cefp1*Ss_zf[ir+1] + cefp2*Ss_zf[ir+2]
    dWEdt_nz[ir] = cefm2*WE_nz[ir-2] + cefm1*WE_nz[ir-1] + cefp0*WE_nz[ir  ] + cefp1*WE_nz[ir+1] + cefp2*WE_nz[ir+2]
    dWEdt_zf[ir] = cefm2*WE_zf[ir-2] + cefm1*WE_zf[ir-1] + cefp0*WE_zf[ir  ] + cefp1*WE_zf[ir+1] + cefp2*WE_zf[ir+2]
  }


  printf( "%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s\n",  \
          "#            time", "dSsdt_nz", "dSsdt_zf", "dWEdt_nz", "dWEdt_zf",  \
                               "RE_nz", "RE_zf", "NE_nz", "NE_zf",  \
                               "Ds_nz", "Ds_zf", "QE" )
  for(ir=1; ir<=2; ir++){
    printf( "%17.7e%17s%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], "NaN", "NaN", "NaN", "NaN",  \
                   RE_nz[ir], RE_zf[ir], NE_nz[ir], NE_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], QE[ir] )
  }
  for(ir=3; ir<=NR-2; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], dSsdt_nz[ir], dSsdt_zf[ir], dWEdt_nz[ir], dWEdt_zf[ir],  \
                   RE_nz[ir], RE_zf[ir], NE_nz[ir], NE_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], QE[ir] )
  }
  for(ir=NR-1; ir<=NR; ir++){
    printf( "%17.7e%17s%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], "NaN", "NaN", "NaN", "NaN",  \
                   RE_nz[ir], RE_zf[ir], NE_nz[ir], NE_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], QE[ir] )
  }
}
