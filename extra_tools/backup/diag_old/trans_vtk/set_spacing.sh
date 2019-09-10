#!/bin/sh

nxw=288
nyw=48
global_nz=32
lx=64.10256410256412
ly=62.83185307179588
lz=3.141592653589793
r_minor=84.0
eps_r=0.18
q_0=1.4
r_major=`echo "scale=14; ${r_minor} / ${eps_r}" | bc`

dx=0`echo "scale=14; ${lx} / ${nxw}" | bc`
dy=`echo "scale=14; ${ly} / ${nyw}" | bc`
dz=`echo "scale=14; ${q_0} * ${r_major} * ${lz} / ${global_nz}" | bc`


loopsta=0
loopend=1300
loopskip=1


echo "Begin phiinxyz."
  loop=${loopsta}
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    input="./stp_phi/phiinxyz_t${cloop}.vtk"
    output=${input}
    sed -e "7c SPACING ${dx} ${dy} ${dz}" ${input} > wk.dat
    mv wk.dat ${output}
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

