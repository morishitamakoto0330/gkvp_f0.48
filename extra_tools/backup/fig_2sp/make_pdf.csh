#!/bin/csh


#### parameter setting ###
set R0_Ln = (2.2 2.2)
set R0_Lt = (6.82 6.82)
set tau = (1.0 1.0)


#### rearrange for plots ###
rm -f ./data/*.dat
./src/calc_elt.sh
cp ../hst/gkvp_f0.48.frq.001 ./data/frq.dat
cat ../hst/gkvp_f0.48.dtc.* > ./data/dtc.dat
cat ../hst/gkvp_f0.48.eng.* > ./data/eng.dat
cat ../hst/gkvp_f0.48.men.* > ./data/men.dat
cat ../hst/gkvp_f0.48.bln.0.* > ./data/bln.0.dat
cat ../hst/gkvp_f0.48.bln.1.* > ./data/bln.1.dat
cat ./src/calc_entropybalance.awk | sed -e "s/R0_Ln =/R0_Ln = $R0_Ln[1]/g" \
                                  | sed -e "s/R0_Lt =/R0_Lt = $R0_Lt[1]/g" \
                                  | sed -e "s/tau =/tau = $tau[1]/g" > wk.0.awk
gawk -f wk.0.awk ./data/bln.0.dat > ./data/ent.0.dat
cat ./src/calc_entropybalance.awk | sed -e "s/R0_Ln =/R0_Ln = $R0_Ln[2]/g" \
                                  | sed -e "s/R0_Lt =/R0_Lt = $R0_Lt[2]/g" \
                                  | sed -e "s/tau =/tau = $tau[2]/g" > wk.1.awk
gawk -f wk.1.awk ./data/bln.1.dat > ./data/ent.1.dat
rm wk.0.awk wk.1.awk


#### plots ###
gnuplot ./src/plot_elt.gn        # Histogram of elapsed time
gnuplot ./src/plot_tdtc.gn       # Time evolution of time step size dt
#gnuplot ./src/plot_grow.gn       # ( For linear runs ) Growthrate for kx = 0 modes
#gnuplot ./src/plot_freq.gn       # ( For linear runs ) Frequency for kx = 0 modes
gnuplot ./src/plot_teng.gn       # Time evolution of electric energy |phi|^2
gnuplot ./src/plot_tmen.gn       # Time evolution of magnetic energy |Al|^2
gnuplot ./src/plot_tent.gn       # Time evolution of entropy balance
gnuplot ./src/plot_tent.0.gn     # Time evolution of entropy balance for electrons
gnuplot ./src/plot_tent.1.gn     # Time evolution of entropy balance for ions
gnuplot ./src/plot_teint.gn      # Time evolution of entropy balance for electric field
gnuplot ./src/plot_tmint.gn      # Time evolution of entropy balance for magnetic field
gnuplot ./src/plot_tent.0.zf.gn  # Time evolution of entropy balance for electrons ( ky = 0 )
gnuplot ./src/plot_tent.1.zf.gn  # Time evolution of entropy balance for ions ( ky = 0 )
gnuplot ./src/plot_teint.zf.gn   # Time evolution of entropy balance for electric field ( ky = 0 )
gnuplot ./src/plot_tmint.zf.gn   # Time evolution of entropy balance for magnetic field ( ky = 0 )
gnuplot ./src/plot_tent.0.nz.gn  # Time evolution of entropy balance for electrons ( ky /= 0 )
gnuplot ./src/plot_tent.1.nz.gn  # Time evolution of entropy balance for ions ( ky /= 0 )
gnuplot ./src/plot_teint.nz.gn   # Time evolution of entropy balance for electric field ( ky /= 0 )
gnuplot ./src/plot_tmint.nz.gn   # Time evolution of entropy balance for magnetic field ( ky /= 0 )


### pdf ###
echo `pwd`'\\\\' | sed -e 's/_/\\_/g' >> parameters.tex
echo `date`'\\\\' >> parameters.tex
cat ../gkvp_f0.48_namelist.001 | sed -e 's/&cmemo/\\underline{cmemo}\\\\/g' \
                               | sed -e 's/&calct/\\underline{calct}\\\\/g' \
                               | sed -e 's/&equib/\\underline{equib}\\\\/g' \
                               | sed -e 's/&run_n/\\underline{run_n}\\\\/g' \
                               | sed -e 's/&files/\\underline{files}\\\\/g' \
                               | sed -e 's/&runlm/\\underline{runlm}\\\\/g' \
                               | sed -e 's/&times/\\underline{times}\\\\/g' \
                               | sed -e 's/&deltt/\\underline{deltt}\\\\/g' \
                               | sed -e 's/&physp/\\underline{physp}\\\\/g' \
                               | sed -e 's/&nperi/\\underline{nperi}\\\\/g' \
                               | sed -e 's/&confp/\\underline{confp}\\\\/g' \
                               | sed -e 's/&vmecp/\\underline{vmecp}\\\\/g' \
                               | sed -e 's/&vmecf/\\underline{vmecf}\\\\/g' \
                               | sed -e 's/&newbz/\\underline{newbz}\\\\/g' \
                               | sed -e 's/&zahyoin/\\underline{zahyoin}\\\\/g' \
                               | sed -e 's/&igsp/\\underline{igsp}\\\\/g' \
                               | sed -e 's/&igsf/\\underline{igsf}\\\\/g' \
                               | sed -e 's/&nu_ref/\\underline{nu_ref}\\\\/g' \
                               | sed -e 's/&end//g' \
                               | sed -e 's/_/\\_/g' \
                               | sed -e 's/$/\\\\/g' >> parameters.tex
echo "#" >> wk.txt
sed -n '/nxw, nyw  =/p'            ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/global_ny =/p'            ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/global_nz =/p'            ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/global_nv, global_nm =/p' ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/nx, ny, nz   =/p'         ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/nv, nm       =/p'         ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/number of species  =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/nproc/p'                  ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
echo "#" >> wk.txt
sed -n '/lx, ly, lz   =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/lz,   z0     =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/lz_l, z0_l   =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/kxmin, kymin =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/kxmax, kymax =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/kperp_max    =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dz           =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dv, vmax     =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dm, mmax     =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
echo "#" >> wk.txt
sed -n '/courant num. =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_perp      =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_zz        =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_vl        =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_col       =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_linear    =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt_max       =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
sed -n '/dt           =/p'   ../log/gkvp_f0.48.000000.0.log.001 >> wk.txt
 
cat wk.txt | sed -e 's/_/\\_/g' | sed -e 's/#/\\#/g' | sed -e 's/$/\\\\/g' >> log.tex
mv parameters.tex pdf/
mv log.tex pdf/
rm wk.txt
cp ./eps/*.eps ./pdf/
cd ./pdf/
latex fig.tex
dvipdf -sPAPERSIZE=a4 fig.dvi
rm fig.aux fig.log fig.dvi
rm *.eps
#gv fig.pdf
