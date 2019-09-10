#!/bin/awk
# <Note>
#   Time derivative is calculated by 4th order central finite difference.
# <How to use>
#   awk -f calc_entropybalance.awk filename > out.dat
#

BEGIN{
  R0_Ln =
  R0_Lt =
  tau =
}

{
      t[NR] = $1                              # Time
  Ss_nz[NR] = $2                              #      Gyrocenter entropy S_s
  Ss_zf[NR] = $3                              #
  WE_nz[NR] = $4                              #         Electric energy W_E
  WE_zf[NR] = $5                              #
  WM_nz[NR] = $6                              #         Magnetic energy W_M
  WM_zf[NR] = $7                              #
  RE_nz[NR] = $8                              #    Electric interaction R_sE
  RE_zf[NR] = $9                              #
  RM_nz[NR] = $10                             #    Magnetic interaction R_sM
  RM_zf[NR] = $11                             #
  NE_nz[NR] = $12                             #   Nonlinear interaction N_sE
  NE_zf[NR] = $13                             #
  NM_nz[NR] = $14                             #   Nonlinear interaction N_sM
  NM_zf[NR] = $15                             #
  Ds_nz[NR] = $16                             # Collisional dissipation D_s
  Ds_zf[NR] = $17                             #
     GE[NR] = tau * ( R0_Ln + R0_Lt ) * $18     #  Electric particle flux T_s * G_sE / L_ps
     GM[NR] = tau * ( R0_Ln + R0_Lt ) * $19     #  Magnetic particle flux T_s * G_sM / L_ps
     QE[NR] = R0_Lt * ( $20 - 2.5 * tau * $18 ) #      Electric heat flux Q_sE
     QM[NR] = R0_Lt * ( $21 - 2.5 * tau * $19 ) #      Magnetic heat flux Q_sM
}

END{
  for(ir=3; ir<=NR-2; ir++){
### for uniform time steps ###
#    dtout_eng =
#    cef = 1 / ( 12 * dtout_eng )
#    dSsdt[ir] = cef * ( - Ss[ir+2] + 8 * Ss[ir+1] - 8 * Ss[ir-1] + Ss[ir-2] )
#    dSedt[ir] = cef * ( - Se[ir+2] + 8 * Se[ir+1] - 8 * Se[ir-1] + Se[ir-2] )
#    dWEdt[ir] = cef * ( - WE[ir+2] + 8 * WE[ir+1] - 8 * WE[ir-1] + WE[ir-2] )
#    dWMdt[ir] = cef * ( - WM[ir+2] + 8 * WM[ir+1] - 8 * WM[ir-1] + WM[ir-2] )
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
    dWMdt_nz[ir] = cefm2*WM_nz[ir-2] + cefm1*WM_nz[ir-1] + cefp0*WM_nz[ir  ] + cefp1*WM_nz[ir+1] + cefp2*WM_nz[ir+2]
    dWMdt_zf[ir] = cefm2*WM_zf[ir-2] + cefm1*WM_zf[ir-1] + cefp0*WM_zf[ir  ] + cefp1*WM_zf[ir+1] + cefp2*WM_zf[ir+2]
  }


  printf( "%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s\n",  \
          "#            time", "dSsdt_nz", "dSsdt_zf", "dWEdt_nz", "dWEdt_zf", "dWMdt_nz", "dWMdt_zf",  \
                               "RE_nz", "RE_zf", "RM_nz", "RM_zf", "NE_nz", "NE_zf", "NM_nz", "NM_zf",  \
                               "Ds_nz", "Ds_zf", "GE", "GM", "QE", "QM" )
  for(ir=1; ir<=2; ir++){
    printf( "%17.7e%17s%17s%17s%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], "NaN", "NaN", "NaN", "NaN", "NaN", "NaN",  \
                   RE_nz[ir], RE_zf[ir], RM_nz[ir], RM_zf[ir], NE_nz[ir], NE_zf[ir], NM_nz[ir], NM_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], GE[ir], GM[ir], QE[ir], QM[ir] )
  }
  for(ir=3; ir<=NR-2; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], dSsdt_nz[ir], dSsdt_zf[ir], dWEdt_nz[ir], dWEdt_zf[ir], dWMdt_nz[ir], dWMdt_zf[ir],  \
                   RE_nz[ir], RE_zf[ir], RM_nz[ir], RM_zf[ir], NE_nz[ir], NE_zf[ir], NM_nz[ir], NM_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], GE[ir], GM[ir], QE[ir], QM[ir] )
  }
  for(ir=NR-1; ir<=NR; ir++){
    printf( "%17.7e%17s%17s%17s%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            t[ir], "NaN", "NaN", "NaN", "NaN", "NaN", "NaN",  \
                   RE_nz[ir], RE_zf[ir], RM_nz[ir], RM_zf[ir], NE_nz[ir], NE_zf[ir], NM_nz[ir], NM_zf[ir],  \
                   Ds_nz[ir], Ds_zf[ir], GE[ir], GM[ir], QE[ir], QM[ir] )
  }
}
