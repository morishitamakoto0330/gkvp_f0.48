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
     Si[NR] = $2      #                 Ion gyrocenter entropy S_i
     WE[NR] = $3      #            Electron gyrocenter entropy S_e
     WM[NR] = $4      #                        Electric energy W_E
   Gies[NR] = $5      #        Electrostatic ion particle flux T_i*(1+eta_i)*G_i
   Giem[NR] = $6      #      Electromagnetic ion particle flux T_i*(1+eta_i)*G_i
   Qies[NR] = $7      #            Electrostatic ion heat flux eta_i*Q_i
   Qiem[NR] = $8      #          Electromagnetic ion heat flux eta_i*Q_i
     Di[NR] = $9      #            Ion collisional dissipation eta_i*Q_i
  Edgei[NR] = $10     #                Ion z-edge compensation Edge_i
    Aiz[NR] = $11     #       Ion artificial dissipation in zz A_i
    Aiv[NR] = $12     #       Ion artificial dissipation in vl A_i
}

END{
  for(ir=3; ir<=NR-2; ir++){
    cef = 1 / ( 12 * dtout_eng )
    dSidt[ir] = cef * ( - Si[ir+2] + 8 * Si[ir+1] - 8 * Si[ir-1] + Si[ir-2] )
    dWEdt[ir] = cef * ( - WE[ir+2] + 8 * WE[ir+1] - 8 * WE[ir-1] + WE[ir-2] )
    dWMdt[ir] = cef * ( - WM[ir+2] + 8 * WM[ir+1] - 8 * WM[ir-1] + WM[ir-2] )
  }

  printf( "%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s%17s\n",  \
          "#            time", "dSidt", "dWEdt", "dWMdt",   \
                               "Gies", "Giem",       \
                               "Qies", "Qiem",       \
                               "Di", "Edgei",        \
                               "Aiz", "Aiv" )
  for(ir=1; ir<=2; ir++){
    printf( "%17.7e%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            time[ir], "NaN", "NaN", "NaN",                 \
                      Gies[ir], Giem[ir],     \
                      Qies[ir], Qiem[ir],     \
                      Di[ir], Edgei[ir],       \
                      Aiz[ir], Aiv[ir] )
  }
  for(ir=3; ir<=NR-2; ir++){
    printf( "%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            time[ir], dSidt[ir], dWEdt[ir], dWMdt[ir], \
                      Gies[ir], Giem[ir],     \
                      Qies[ir], Qiem[ir],     \
                      Di[ir], Edgei[ir],       \
                      Aiz[ir], Aiv[ir] )
  }
  for(ir=NR-1; ir<=NR; ir++){
    printf( "%17.7e%17s%17s%17s%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e%17.7e\n",  \
            time[ir], "NaN", "NaN", "NaN",                 \
                      Gies[ir], Giem[ir],     \
                      Qies[ir], Qiem[ir],     \
                      Di[ir], Edgei[ir],       \
                      Aiz[ir], Aiv[ir] )
  }
}
