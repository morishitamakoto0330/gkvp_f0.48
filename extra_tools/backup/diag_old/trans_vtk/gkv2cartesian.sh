#!/bin/sh

r_minor=84.0     # Minor radius [in rho unit]
eps_r=0.18       # Inverse aspect ratio (=r/R_0)
q_0=1.4          # Safety factor
s_hat=0.78       # Magnetic shear


loopsta=0
loopend=117
loopskip=1


echo "Begin phiinxyz."
  loop=${loopsta}
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    input="./stp_gkv/phiinxyz_t${cloop}.vtk"
    output="./stg_cartesian/phiinxyz_t${cloop}.vtk"
    awk -f ./gkv2cartesian.awk \
        -v r_minor=${r_minor} eps_r=${eps_r}  \
           q_0=${q_0} s_hat=${s_hat} \
        ${input} > ${output}
    loop=`expr ${loop} + ${loopskip}`
  done

#echo "Begin Alinxyz."
#  loop=${loopsta}
#  while [ $loop -le ${loopend} ]
#  do
#    cloop=`printf "%08d" $loop`
#    input="./stp_gkv/Alinxyz_t${cloop}.vtk"
#    output="./stg_cartesian/Alinxyz_t${cloop}.vtk"
#    awk -f ./gkv2cartesian.awk \
#        -v r_minor=${r_minor} eps_r=${eps_r}  \
#           q_0=${q_0} s_hat=${s_hat} \
#        ${input} > ${output}
#    loop=`expr ${loop} + ${loopskip}`
#  done
#
#
#  is=0
#  while [ $is -le 1 ]
#  do
#    cis=`printf "%01d" $is`
#
#echo "Begin densinxyz."
#    loop=${loopsta}
#    while [ $loop -le ${loopend} ]
#    do
#      cloop=`printf "%08d" $loop`
#      input="./stp_gkv/densinxyz_ns${cis}_t${cloop}.vtk"
#      output="./stg_cartesian/densinxyz_ns${cis}_t${cloop}.vtk"
#      awk -f ./gkv2cartesian.awk \
#          -v r_minor=${r_minor} eps_r=${eps_r}  \
#             q_0=${q_0} s_hat=${s_hat} \
#          ${input} > ${output}
#    loop=`expr ${loop} + ${loopskip}`
#    done
#  
#echo "Begin uparainxyz."
#    loop=${loopsta}
#    while [ $loop -le ${loopend} ]
#    do
#      cloop=`printf "%08d" $loop`
#      input="./stp_gkv/uparainxyz_ns${cis}_t${cloop}.vtk"
#      output="./stg_cartesian/uparainxyz_ns${cis}_t${cloop}.vtk"
#      awk -f ./gkv2cartesian.awk \
#          -v r_minor=${r_minor} eps_r=${eps_r}  \
#             q_0=${q_0} s_hat=${s_hat} \
#          ${input} > ${output}
#    loop=`expr ${loop} + ${loopskip}`
#    done
#  
#echo "Begin pressianxyz."
#    loop=${loopsta}
#    while [ $loop -le ${loopend} ]
#    do
#      cloop=`printf "%08d" $loop`
#      input="./stp_gkv/presinxyz_ns${cis}_t${cloop}.vtk"
#      output="./stg_cartesian/presinxyz_ns${cis}_t${cloop}.vtk"
#      awk -f ./gkv2cartesian.awk \
#          -v r_minor=${r_minor} eps_r=${eps_r}  \
#             q_0=${q_0} s_hat=${s_hat} \
#          ${input} > ${output}
#    loop=`expr ${loop} + ${loopskip}`
#    done
#  
#echo "Begin qparaianxyz."
#    loop=${loopsta}
#    while [ $loop -le ${loopend} ]
#    do
#      cloop=`printf "%08d" $loop`
#      input="./stp_gkv/qparainxyz_ns${cis}_t${cloop}.vtk"
#      output="./stg_cartesian/qparainxyz_ns${cis}_t${cloop}.vtk"
#      awk -f ./gkv2cartesian.awk \
#          -v r_minor=${r_minor} eps_r=${eps_r}  \
#             q_0=${q_0} s_hat=${s_hat} \
#          ${input} > ${output}
#    loop=`expr ${loop} + ${loopskip}`
#    done
#
#    is=`expr ${is} + 1`
#  done

