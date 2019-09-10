#!/bin/sh

r_minor=84.0     # Minor radius [in rho unit]
eps_r=0.18       # Inverse aspect ratio (=r/R_0)
q_0=1.4          # Safety factor
s_hat=0.78       # Magnetic shear
n_alp=12

n1=41
n2=41
n3=65



loopsta=0
loopend=117
loopskip=1


echo "Begin phiinxyz."
  loop=${loopsta}
  while [ $loop -le ${loopend} ]
  do
    cloop=`printf "%08d" $loop`
    input="./stp_gkv/phiinxyz_t${cloop}.vtk"
    output="./stg_full/phiinxyz_inside_t${cloop}.vtk"

    awk -f ./gkv2full1_surf.awk \
        -v r_minor=${r_minor} eps_r=${eps_r}  \
           q_0=${q_0} s_hat=${s_hat} n_alp=${n_alp} \
        ${input} > ${output}

    i3=0
    while [ $i3 -lt $n3 ]
    do
      awk -f ./gkv2full2_surf.awk \
          -v n_alp=${n_alp} n1=${n1} n2=${n2} i3=${i3} \
          ${input} >> ${output}
      i3=`expr $i3 + 1`
    done

    loop=`expr ${loop} + ${loopskip}`
  done
