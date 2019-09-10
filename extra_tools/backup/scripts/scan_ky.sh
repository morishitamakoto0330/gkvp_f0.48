#!/bin/sh

#for param in `seq 20 10 100`; do # ky scan  (ky = param1 * 1d-2; n_alp = param1; q_d = 1/(2*param1))
for param in `seq 100 100 2000`; do # ky scan  (ky = param1 * 1d-2; n_alp = param1; q_d = 1/(2*param1))

  CURRENTDIR=`pwd`
 
  ### copy run dir ###                                                                                 ### for debug ###
  OLD_RUNDIR="ky0010"
  NEW_RUNDIR=`echo ${param} | awk '{printf("%s%04d","ky",$1)}'`;                                       echo ${NEW_RUNDIR}
  mkdir -p ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/gkvp_mpifft.exe         ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/go.sh                   ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/sub.q                   ${NEW_RUNDIR}/run
  cp ${OLD_RUNDIR}/run/gkvp_f0.48_namelist.001 ${NEW_RUNDIR}/run
 
  ### change work dir name ###
  OLD_WKDIR="\/large2\/z41049a\/gkvp\/f0.48\/applegate2007ppcf_salpha_lb\/linear_kdep\/ky0010"
  NEW_WKDIR="\/large2\/z41049a\/gkvp\/f0.48\/applegate2007ppcf_salpha_lb\/linear_kdep\/${NEW_RUNDIR}"; #echo ${NEW_WKDIR}
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/go.sh
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/sub.q
  sed -i -e "s/${OLD_WKDIR}/${NEW_WKDIR}/g" ${NEW_RUNDIR}/run/gkvp_f0.48_namelist.001
 
  ### change parameter ###
  n_alp=`echo ${param} | awk '{printf("%d",$1)}'`
  OLD_PARAM1="n_alp=10"
  NEW_PARAM1="n_alp=${n_alp}";                                                                         #echo ${NEW_PARAM1}
  q_d=`echo ${param} | awk '{printf("%15.10f",1.0/(2.0*$1*1))}'`
  OLD_PARAM2="q_d      = 0.05d0"
  NEW_PARAM2="q_d      = ${q_d}d0";                                                                    #echo ${NEW_PARAM2}
  sed -i -e "s/${OLD_PARAM1}/${NEW_PARAM1}/g" ${NEW_RUNDIR}/run/gkvp_f0.48_namelist.001
  sed -i -e "s/${OLD_PARAM2}/${NEW_PARAM2}/g" ${NEW_RUNDIR}/run/gkvp_f0.48_namelist.001
 
  ### execution ###                                                                                    #################
  cd ${NEW_RUNDIR}/run
  #./go.sh
  
  cd ${CURRENTDIR}

done
