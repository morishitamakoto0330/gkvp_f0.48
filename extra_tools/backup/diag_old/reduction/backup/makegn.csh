#!/bin/csh

### DIAG parameters ###
  set time_spacing = 10

  set flag_cnt = 0
  set flag_ffinzv      = 0
  set flag_mom = 1
  set flag_phiinkxky   = 1
  set flag_phiinxy     = 1
  set flag_phiinz      = 1
  set flag_Alinkxky    = 1
  set flag_Alinxy      = 1
  set flag_Alinz       = 1
  set flag_fluxinkxky  = 1
  set flag_fluxinky    = 1
  set flag_transinkxky = 1
#######################

### pm3d palette sample ###
#  echo "  set palette rgbformulae 23,22,21  # cold"
#  echo "  set palette rgbformulae 21,22,23  # hot"
#  echo "  set palette rgbformulae 22,13,-31 # rainbow"
#  echo "  set palette define ( -1 'blue', 0 'white', 1 'red' ) # positive and negative"
###########################


  if ( $flag_cnt ) then

    if ( $flag_ffinzv ) then
      set plotfile = plot_ffinzv.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set contour base"      >> $plotfile
      echo "  set cntrparam levels 20"   >> $plotfile
      echo "  unset surface"         >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set view 0, 0"         >> $plotfile
      echo "  set xlabel 'Filed-aligned coordinate z'"   >> $plotfile
      echo "  set ylabel 'Parallel velocity vl'" >> $plotfile
      awk -f makegn_xytime.awk ./plt/ffinzvtime.dat time_spacing=${time_spacing} >> $plotfile
    endif

  endif

  if ( $flag_mom ) then

    if ( $flag_phiinxy ) then
      set plotfile = plot_phiinxy.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial direction x'"   >> $plotfile
      echo "  set ylabel 'Poloidal direction y'" >> $plotfile
      echo "  set palette rgbformulae 23,22,21  # cold" >> $plotfile
      awk -f makegn_xytime.awk ./plt/phiinxytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_phiinkxky ) then
      set plotfile = plot_phiinkxky.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial wave number kx'"   >> $plotfile
      echo "  set ylabel 'Poloidal wave number ky'" >> $plotfile
      echo "  set palette rgbformulae 23,22,21  # cold" >> $plotfile
      awk -f makegn_xytime.awk ./plt/phiinkxkytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_phiinz ) then
      set plotfile = plot_phiinz.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Field-aligned coordinate z'"           >> $plotfile
      echo "  set ylabel 'Re[{/Symbol f}_k], Im[{/Symbol f}_k]'" >> $plotfile
      awk -f makegn_ztime.awk ./plt/phiinztime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_Alinxy ) then
      set plotfile = plot_Alinxy.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial direction x'"   >> $plotfile
      echo "  set ylabel 'Poloidal direction y'" >> $plotfile
      echo "  set palette rgbformulae 21,22,23  # hot" >> $plotfile
      awk -f makegn_xytime.awk ./plt/Alinxytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_Alinkxky ) then
      set plotfile = plot_Alinkxky.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial wave number kx'"   >> $plotfile
      echo "  set ylabel 'Poloidal wave number ky'" >> $plotfile
      echo "  set palette rgbformulae 21,22,23  # hot" >> $plotfile
      awk -f makegn_xytime.awk ./plt/Alinkxkytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_Alinz ) then
      set plotfile = plot_Alinz.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Field-aligned coordinate z'"           >> $plotfile
      echo "  set ylabel 'Re[{/Symbol f}_k], Im[{/Symbol f}_k]'" >> $plotfile
      awk -f makegn_ztime.awk ./plt/Alinztime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_fluxinkxky ) then
      set plotfile = plot_fluxinkxky.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial wave number kx'"   >> $plotfile
      echo "  set ylabel 'Poloidal wave number ky'" >> $plotfile
      echo "  set palette rgbformulae 22,13,-31  # rainbow" >> $plotfile
      awk -f makegn_xytime.awk ./plt/fluxinkxkytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_fluxinky ) then
      set plotfile = plot_fluxinkytime.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Poloidal wave number ky'" >> $plotfile
      echo "  set ylabel 'Time t'"   >> $plotfile
      echo "  set palette rgbformulae 22,13,-31  # rainbow" >> $plotfile
      echo "  splot './plt/fluxinkytime.dat' u 1:2:3" >> $plotfile
      echo "  pause -1" >> $plotfile
      set plotfile = plot_fluxinky.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Poloidal wave number ky'"           >> $plotfile
      echo "  set ylabel 'Flux'" >> $plotfile
      awk -f makegn_ztime.awk ./plt/fluxinkytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

    if ( $flag_transinkxky ) then
      set plotfile = plot_transinkxky.gn
      echo "#\!/usr/bin/gnuplot"  > $plotfile
      echo "# set terminal jpeg" >> $plotfile
      echo "  set pm3d map"      >> $plotfile
      echo "  set size square"   >> $plotfile
      echo "  unset key"         >> $plotfile
      echo "  set xlabel 'Radial wave number kx'"   >> $plotfile
      echo "  set ylabel 'Poloidal wave number ky'" >> $plotfile
      echo "  set palette define ( -1 'blue', 0 'white', 1 'red' ) # positive and negative" >> $plotfile
      awk -f makegn_xytime.awk ./plt/transinkxkytime.dat time_spacing=${time_spacing} >> $plotfile
    endif

  endif
  
