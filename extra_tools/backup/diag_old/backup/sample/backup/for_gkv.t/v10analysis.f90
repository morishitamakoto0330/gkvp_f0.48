MODULE analysis_header
!-------------------------------------------------------------------------------
!
!     Header for data analysis
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  implicit none


  integer, parameter :: snum = 1      ! begining of simulation runs
  integer, parameter :: enum = 1      ! end of simulation runs


  integer, parameter :: DP = selected_real_kind(14)

  integer :: ns, nxw, nyw, global_nz, global_nv, global_nm,           &
             nprocz, nprocv, nprocm, nx, ny, nz, nv, nm,              &
             n_alp, n_tht

  real(kind=DP) :: eta_i, eta_e, mei, tau, beta, lambda2, nu_i, nu_e, &
                   diff_z, diff_v,                                    &
                   rad_a, eps_b, eps_r, r_minor, r_major, rho2L_n,    &
                   q_0, q_d, s_hat,                                   &
                   eps_h, lprd, mprd, malpha,                         &
                   lx, ly, lz, vmax, mmax,                            &
                   kxmin, kymin, dz, dpara, dv, dm

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, twopi = pi * 2._DP
  real(kind=DP),    parameter :: eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )

  real(kind=DP), dimension(:), allocatable :: m_s, e_s, T_s, eta_s
  real(kind=DP), dimension(:), allocatable :: xx, yy, kx, ky, gzz, gvl, gmu

  complex(kind=DP), dimension(:), allocatable :: ck
  integer, dimension(:),          allocatable :: dj
  real(kind=DP), dimension(:), allocatable :: gbeq
  real(kind=DP), dimension(:,:,:), allocatable :: gfmx
  real(kind=DP) :: cfsrf

! --- unit numbers for I/O
  integer, parameter :: olog = 10, &
                        iset = 12, &
                        ofzv = 20, pfzv = 21, &
                        ozzv = 30, pzzv = 31, &
                        ozvm = 32, pzvm = 33, &
                        opkk = 40, ppkk = 41, &
                        opxy = 42, ppxy = 43, &
                        opzz = 44, ppzz = 45, &
                        oakk = 50, pakk = 51, &
                        oaxy = 52, paxy = 53, &
                        oazz = 54, pazz = 55, &
                        oxtk = 60, pxtk = 61, &
                        oxkk = 62, pxkk = 63, &
                        oxtx = 64, pxtx = 65

 CONTAINS

SUBROUTINE set_parameters( inum )

  implicit none

  integer, intent(in) :: inum

  character*3 :: cnum

    write( cnum, fmt="(i3.3)" ) inum

    open( iset, file="./bdata/gkv.set."//cnum, status="old",  &
          action="read", form="unformatted" )

    read( iset ) ns, nxw, nyw, global_nz, global_nv, global_nm,     &
                 nprocz, nprocv, nprocm, nx, ny, nz, nv, nm,        &
                 eta_i, eta_e, mei, tau, beta, lambda2, nu_i, nu_e, &
                 diff_z, diff_v,                                    &
                 rad_a, eps_b, eps_r, r_minor, r_major, rho2L_n,    &
                 q_0, q_d, s_hat,                                   &
                 eps_h, lprd, mprd, malpha,                         &
                 n_alp, n_tht, lx, ly, lz, vmax, mmax,              &
                 kxmin, kymin, dz, dpara, dv, dm
    close( iset )

    allocate( m_s(1:ns) )
    allocate( e_s(1:ns) )
    allocate( T_s(1:ns) )
    allocate( eta_s(1:ns) )
    allocate( xx(0:2*nxw-1) )
    allocate( yy(0:2*nyw-1) )
    allocate( kx(0:nx) )
    allocate( ky(-ny:ny) )
    allocate( gzz(-global_nz:global_nz-1) )
    allocate( gvl(1:2*global_nv) )
    allocate( gmu(0:global_nm) )
    allocate( ck(-ny:ny) )
    allocate( dj(-ny:ny) )
    allocate( gbeq(-global_nz:global_nz-1) )
    allocate( gfmx(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) )


END SUBROUTINE set_parameters

SUBROUTINE reset_parameters

  implicit none

    deallocate( m_s )
    deallocate( e_s )
    deallocate( T_s )
    deallocate( eta_s )
    deallocate( xx )
    deallocate( yy )
    deallocate( kx )
    deallocate( ky )
    deallocate( gzz )
    deallocate( gvl )
    deallocate( gmu )
    deallocate( ck )
    deallocate( dj )
    deallocate( gbeq )
    deallocate( gfmx )


END SUBROUTINE reset_parameters

END MODULE analysis_header


MODULE analysis_fft
!-------------------------------------------------------------------------------
!
!    FFT module using fftw3
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

  include "fftw3.f"

  private

  integer(kind=DP), save :: plan

  public   fft_pre, fft_backward


CONTAINS


!--------------------------------------
  SUBROUTINE fft_pre
!--------------------------------------
!  Initialization of FFT

    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wwxy


      call dfftw_plan_dft_c2r_2d( plan, 2*nxw, 2*nyw, &
                                  wwkk, wwxy,       &
                                  FFTW_ESTIMATE )


  END SUBROUTINE fft_pre


!--------------------------------------
  SUBROUTINE fft_backward ( gww, wwxy )
!--------------------------------------
!  Execution of FFT

    complex(kind=DP), intent(in), &
      dimension(0:nx,-ny:ny)                     :: gww
    real(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:2*nyw-1)             :: wwxy

    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
    integer :: mx, my, iz


      do my = 0, 2*nyw-1
        do mx = 0, nxw
          wwkk(mx,my) = ( 0._DP, 0._DP )
        end do
      end do
      do my = 0, ny
        do mx = 0, nx
          wwkk(mx,my) = gww(mx,my)
        end do
      end do
      do my = -ny, -1
        do mx = 0, nx
          wwkk(mx,2*nyw+my) = gww(mx,my)
        end do
      end do
      mx = 0
        do my = 1, ny
          wwkk(mx,2*nyw-my) = conjg( wwkk(mx,my) )
        end do

    call dfftw_execute_dft( plan, wwkk, wwxy )


  END SUBROUTINE fft_backward


END MODULE analysis_fft


PROGRAM analysis
!-------------------------------------------------------------------------------
!
!     Data analysis from binary output
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header, only : snum, enum, set_parameters, reset_parameters, olog
  use analysis_fft, only : fft_pre

  implicit none

  integer :: inum
  character*3 :: cnum

    call file_open

    do inum = snum, enum

      if ( inum == snum ) then
        call set_parameters( inum )
        call fft_pre
        call init
      end if

      write( cnum, fmt="(i3.3)" ) inum
      write( olog, * ) " ##### inum = "//cnum//" #####"

!      call analysis_cnt( inum )
      call analysis_zvm( inum )
      call analysis_mom( inum )

      if ( inum == enum ) then
        call reset_parameters
      end if

      write( olog, * ) " ######################"
      write( olog, * ) ""

    end do

    call file_close

  stop

  CONTAINS


SUBROUTINE file_open
!-------------------------------------------------------------------------------
!
!     Open files
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

    open( olog, file="./plt/log.dat" )

    open( pfzv, file="./plot_ffinzv.gn" )
      write( pfzv, * ) "set contour base"
      write( pfzv, * ) "set cntrparam levels 20"
      write( pfzv, * ) "set nosurface"
      write( pfzv, * ) "set view 0, 0"
      write( pfzv, * ) "set xlabel 'Field-aligned coordinate z'"
      write( pfzv, * ) "set ylabel 'Parallel velocity vl'"

    open( pzzv, file="./plot_fzvminzv.gn" )
      write( pzzv, * ) "set contour base"
      write( pzzv, * ) "set cntrparam levels 20"
      write( pzzv, * ) "set nosurface"
      write( pzzv, * ) "set view 0, 0"
      write( pzzv, * ) "set xlabel 'Field-aligned coordinate z'"
      write( pzzv, * ) "set ylabel 'Parallel velocity vl'"

    open( pzvm, file="./plot_fzvminvm.gn" )
      write( pzvm, * ) "set contour base"
      write( pzvm, * ) "set cntrparam levels 20"
      write( pzvm, * ) "set nosurface"
      write( pzvm, * ) "set view 0, 0"
      write( pzvm, * ) "set xlabel 'Parallel velocity vl'"
      write( pzvm, * ) "set ylabel 'Magnetic moment mu'"

    open( ppkk, file="./plot_phiinkxky.gn" )
      write( ppkk, * ) "set pm3d map"
      write( ppkk, * ) "set size square"
      write( ppkk, * ) "set xlabel 'Wave number kx'"
      write( ppkk, * ) "set ylabel 'Wave number ky'"

    open( ppxy, file="./plot_phiinxy.gn" )
      write( ppxy, * ) "set pm3d map"
      write( ppxy, * ) "set size square"
      write( ppxy, * ) "set xlabel 'Radial direction x'"
      write( ppxy, * ) "set ylabel 'Poloidal direction y'"

    open( ppzz, file="./plot_phiinzz.gn" )
      write( ppzz, * ) "set xlabel 'Field-aligned coordinate z'"
      write( ppzz, * ) "set ylabel '| phi |'"

    open( pakk, file="./plot_Alinkxky.gn" )
      write( pakk, * ) "set pm3d map"
      write( pakk, * ) "set size square"
      write( pakk, * ) "set xlabel 'Wave number kx'"
      write( pakk, * ) "set ylabel 'Wave number ky'"

    open( paxy, file="./plot_Alinxy.gn" )
      write( paxy, * ) "set pm3d map"
      write( paxy, * ) "set size square"
      write( paxy, * ) "set xlabel 'Radial direction x'"
      write( paxy, * ) "set ylabel 'Poloidal direction y'"

    open( pazz, file="./plot_Alinzz.gn" )
      write( pazz, * ) "set xlabel 'Field-aligned coordinate z'"
      write( pazz, * ) "set ylabel '| Al |'"

    open( oxtk, file="./plt/tfluxes.dat" )
      write( oxtk, "(99a17)" ) "#            Time",                   &
                               "p_flux_es(ns)", "p_flux_em(ns)",  &
                               "h_flux_es(ns)", "h_flux_em(ns)"

    open( pxtk, file="./plot_tfluxes.gn" )
      write( pxtk, * ) "set xlabel 'Time t v_{ti}/L_n'"
      write( pxtk, * ) "set ylabel 'Particle Fluxes'"
      write( pxtk, * ) "plot './plt/tfluxes.dat' u 1:2 ti '{/Symbol G}_{iE}' w l, \"
      write( pxtk, * ) "                      '' u 1:3 ti '{/Symbol G}_{eE}' w l, \"
      write( pxtk, * ) "                      '' u 1:4 ti '{/Symbol G}_{iM}' w l, \"
      write( pxtk, * ) "                      '' u 1:5 ti '{/Symbol G}_{eM}' w l"
      write( pxtk, * ) "pause -1"
      write( pxtk, * ) "set ylabel 'Heat Fluxes'"
      write( pxtk, * ) "plot './plt/tfluxes.dat' u 1:6 ti 'Q_{iE}' w l, \"
      write( pxtk, * ) "                      '' u 1:7 ti 'Q_{eE}' w l, \"
      write( pxtk, * ) "                      '' u 1:8 ti 'Q_{iM}' w l, \"
      write( pxtk, * ) "                      '' u 1:9 ti 'Q_{eM}' w l"
      write( pxtk, * ) "pause -1"

    open( pxkk, file="./plot_fluxinkxky.gn" )
      write( pxkk, * ) "set pm3d map"
      write( pxkk, * ) "set size square"
      write( pxkk, * ) "set xlabel 'Wave number kx'"
      write( pxkk, * ) "set ylabel 'Wave number ky'"

    open( oxtx, file="./plt/fluxintx.dat" )
      write( oxtx, "(99a17)" ) "#            Time", "x",          &
                               "p_flux_es(ns)", "p_flux_em(ns)",  &
                               "h_flux_es(ns)", "h_flux_em(ns)"

    open( pxtx, file="./plot_fluxintx.gn" )
      write( pxtx, * ) "set pm3d map"
      write( pxtx, * ) "set xlabel 'Time t v_{ti}/L_n'"
      write( pxtx, * ) "set ylabel 'Radial direction x'"
      write( pxtx, * ) "splot './plt/fluxintx.dat' u 1:2:3"
      write( pxtx, * ) "pause -1"


END SUBROUTINE file_open


SUBROUTINE file_close
!-------------------------------------------------------------------------------
!
!     Close files
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

    close( olog )
    close( pfzv )
    close( pzzv )
    close( pzvm )
    close( ppkk )
    close( ppxy )
    close( ppzz )
    close( pakk )
    close( paxy )
    close( pazz )
    close( oxtk )
    close( pxtk )
    close( pxkk )
    close( oxtx )
    close( pxtx )


END SUBROUTINE file_close


SUBROUTINE init
!-------------------------------------------------------------------------------
!
!     Initial setting
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

  real(kind=DP) :: lmmq
  integer :: mx, my, iz, iv, im, is

!!!! for kinetic ions with adiabatic electron response
!    is = 1
!    m_s(is) = 1._DP
!    e_s(is) = 1._DP
!    T_s(is) = 1._DP
!    eta_s(is) = eta_i

!!!! for kinetic electrons without ion perturbations
!    is = 1
!    m_s(is) = mei
!    e_s(is) = -1._DP
!    T_s(is) = 1._DP / tau
!    eta_s(is) = eta_e

!!!! for kinetic electrons and ions
    is = 1                ! for ion
    m_s(is) = 1._DP
    e_s(is) = 1._DP
    T_s(is) = 1._DP
    eta_s(is) = eta_i
    is = 2                ! for electron
    m_s(is) = mei
    e_s(is) = -1._DP
    T_s(is) = 1._DP / tau
    eta_s(is) = eta_e

    do mx = 0, 2*nxw-1
      xx(mx) = -lx + lx/real(nxw)*mx
    end do

    do my = 0, 2*nyw-1
      yy(my) = -ly + ly/real(nyw)*my
    end do

    do mx = 0, nx
      kx(mx) = kxmin * real( mx, kind=DP )
    end do

    do my = -ny, ny
      ky(my) = kymin * real( my, kind=DP )
    end do

    do iz = -global_nz, global_nz-1
      gzz(iz) = dz * real( iz, kind=DP )
    end do

    do iv = 1, 2*global_nv
      gvl(iv) = - vmax + dv * real( iv-1, kind=DP )
    end do

    do im = 0, global_nm
      gmu(im) = 0.5_DP * ( dm * real( im, kind=DP ) )**2
    end do

    do my = -ny, ny
      ck(my)   = exp( ui * 2._DP * pi * q_0 &
                         * real( n_alp * n_tht * my, kind=DP ) )
      dj(my)   = - nint( real( 2 * n_alp * n_tht * my, kind=DP ) * q_d )
    end do

    lmmq = lprd - mprd * q_0
    do iz = -global_nz, global_nz-1
    !%%% original GKV %%%
    !  gbeq(iz) = 1._DP - eps_r * cos( gzz(iz) ) &
    !                   - eps_h * cos( lmmq * gzz(iz) - malpha )
    !%%% Merz's thesis and GENE source code %%%
      gbeq(iz) = 1._DP / ( 1._DP + eps_r * cos( gzz(iz) ) )
    end do
    cfsrf = 0._DP
    do iz = -global_nz, global_nz-1
      cfsrf = cfsrf + 1._DP / gbeq(iz)
    end do

    do im = 0, global_nm
      do iv = 1, 2*global_nv
        do iz = -global_nz, global_nz-1
          gfmx(iz,iv,im) = exp( - 0.5_DP * gvl(iv)**2 - gbeq(iz) * gmu(im) ) &
                         / sqrt( twopi**3 )
        end do
      end do
    end do


      write( olog, * ) " ##### parameters #####"
      write( olog, * ) ""
      write( olog, * ) " # discretizations "
      write( olog, * ) " # ns           = ", ns
      write( olog, * ) " # nxw          = ", nxw
      write( olog, * ) " # nyw          = ", nyw
      write( olog, * ) " # global_nz    = ", global_nz
      write( olog, * ) " # global_nv    = ", global_nv
      write( olog, * ) " # global_nm    = ", global_nm
      write( olog, * ) " # nprocz       = ", nprocz
      write( olog, * ) " # nprocv       = ", nprocv
      write( olog, * ) " # nprocm       = ", nprocm
      write( olog, * ) " # nx, ny, nz   = ", nx, ny, nz
      write( olog, * ) " # nv, nm       = ", nv, nm
      write( olog, * ) ""
      write( olog, * ) " # physics "
      write( olog, * ) " # eta_i        = ", eta_i
      write( olog, * ) " # eta_e        = ", eta_e
      write( olog, * ) " # mei          = ", mei
      write( olog, * ) " # tau          = ", tau
      write( olog, * ) " # beta         = ", beta
      write( olog, * ) " # lambda2      = ", lambda2
      write( olog, * ) " # nu_i         = ", nu_i
      write( olog, * ) " # nu_e         = ", nu_e
      write( olog, * ) " # diff_z       = ", diff_z
      write( olog, * ) " # diff_v       = ", diff_v
      write( olog, * ) " # rad_a        = ", rad_a
      write( olog, * ) " # eps_b        = ", eps_b
      write( olog, * ) " # eps_r        = ", eps_r
      write( olog, * ) " # r_minor/rho  = ", r_minor
      write( olog, * ) " # r_major/L_n  = ", r_major
      write( olog, * ) " # rho    /L_n  = ", rho2L_n
      write( olog, * ) " # q_0, q_d     = ", q_0, q_d
      write( olog, * ) " # s_hat        = ", s_hat
      write( olog, * ) " # (For a helical system)"
      write( olog, * ) " # eps_h        = ", eps_h
      write( olog, * ) " # lprd         = ", lprd
      write( olog, * ) " # mprd         = ", mprd
      write( olog, * ) " # malpha       = ", malpha
      write( olog, * ) ""
      write( olog, * ) " # domains "
      write( olog, * ) " # n_alp, n_tht = ", n_alp, n_tht
      write( olog, * ) " # lx, ly, lz   = ", lx, ly, lz
      write( olog, * ) " # kxmin, kymin = ", kxmin, kymin
      write( olog, * ) " # dz, dpara    = ", dz, dpara
      write( olog, * ) " # dv, vmax     = ", dv, vmax
      write( olog, * ) " # dm, mmax     = ", dm, mmax
      write( olog, * ) ""
      write( olog, * ) ""


END SUBROUTINE init


SUBROUTINE analysis_cnt( inum )
!-------------------------------------------------------------------------------
!
!     Data analysis from gkv.****.cnt.%%%
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer, intent(in) :: inum

! --- local variables

  character*3 :: cnum
  integer :: ir, irz, irv, irm   ! rank, rankz, rankv, rankm
  character*4 :: crank
  integer :: loop, ios

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-nz:nz-1,1:2*nv,0:nm,1:ns)         :: cf
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm,1:ns) :: gff
  integer :: mx, my, iz, iv, im, is, giz, giv, gim

    write( cnum, fmt="(i3.3)" ) inum
    write( olog, * ) " # Analysis of cnt."

! --- file open ---
    do ir = 0, nprocz * nprocv * nprocm - 1
      write( crank, fmt="(i4.4)" ) ir
      open( unit=1000+ir, file="./bdata/gkv."//crank//".cnt."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".cnt."//cnum, 1000+ir
    end do


! --- time-step loop ---
    loop = 0
    do

  if ( loop >= 0 ) then            !%%% to restart

      do ir = 0, nprocz * nprocv * nprocm - 1

        read( unit=1000+ir, iostat=ios ) time, cf
        if ( ios /= 0 ) exit

        irz = mod( ir, nprocz )
        irv = mod( ir / nprocz, nprocv )
        irm = int( ir / ( nprocz * nprocv ) )
        do is = 1, ns
          do im = 0, nm
            gim = (nm+1) * irm + im
            do iv = 1, 2*nv
              giv = 2*nv * irv + iv
              do iz = -nz, nz-1
                giz = - global_nz + 2*nz * irz + iz + nz
                do my = -ny, ny
                  do mx = 0, nx
                    gff(mx,my,giz,giv,gim,is) = cf(mx,my,iz,iv,im,is)
                  end do
                end do
              end do
            end do
          end do
        end do


      end do

  else                             !%%% to restart

      do ir = 0, nprocz * nprocv * nprocm - 1

        read( unit=1000+ir, iostat=ios )
        if ( ios /= 0 ) exit

      end do

  end if                           !%%% to restart

      if ( ios < 0 ) then
        write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                         ios, 1000+ir
        write( olog, * ) ""
        exit
      else if ( ios > 0 ) then
        write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                         ios, 1000+ir
        write( olog, * ) ""
        exit
      end if

      call ffinzv( inum, loop, time, gff )
      write( olog, * ) " # Data output at loop, time = ", loop, time

      loop = loop + 1

    end do


! --- file close ---
    do ir = 0, nprocz * nprocv * nprocm - 1
      close( unit=1000+ir )
    end do


END SUBROUTINE analysis_cnt


SUBROUTINE ffinzv( inum, loop, time, gff )
!-------------------------------------------------------------------------------
!
!     Write ff in zz, vl space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                       :: inum
  integer,          intent(in)                       :: loop
  real(kind=DP),    intent(in)                       :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm,1:ns) :: gff

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: mx, my, iz, iv, im

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pfzv, * ) "# set output './plt/ffinzv_t"//cloop//"."//cnum//".jpg'"
    write( pfzv, * ) "splot './plt/ffinzv_t"//cloop//"."//cnum//"' u 1:2:3 w pm3d"
    write( pfzv, * ) "pause 0.2"

  if ( loop >= 0 ) then            !%%% to restart

    open( ofzv, file="./plt/ffinzv_t"//cloop//"."//cnum )
    write( ofzv, * ) "#  time = ", time 
    write( ofzv, "(99a17)" ) "#              zz","vl","|ff|"
    mx = 0
    my = 4
    im = 4

    if ( dj(my) == 0 ) then

      do iv = 1, 2*global_nv
        do iz = -global_nz, global_nz-1
          write( ofzv, "(99G17.7E3)" ) gzz(iz), gvl(iv),  &
                                     abs( gff(mx,my,iz,iv,im,1:ns) )
        end do
        write( ofzv, * )
      end do

    else

    do iv = 1, 2*global_nv

      connect_min = int( ( nx + mx ) / abs( dj(my) ) )
      if ( connect_min .ne. 0 ) then
        do iconnect = connect_min, 1, -1
          mxw = mx+iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( ofzv, "(99G17.7E3)" )                               &
                  - twopi * real(iconnect) + gzz(iz), gvl(iv),           &
                  abs( ck(my)**iconnect * conjg( gff(-mxw,-my,iz,iv,im,1:ns) ) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( ofzv, "(99G17.7E3)" )                               &
                  - twopi * real(iconnect) + gzz(iz), gvl(iv),           &
                  abs( ck(my)**iconnect * gff(mxw,my,iz,iv,im,1:ns) )
              end do

          end if

        end do
      end if

      connect_max = int( ( nx - mx ) / abs( dj(my) ) )
        do iconnect = 0, connect_max
          mxw = mx-iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( ofzv, "(99G17.7E3)" )                               &
                  + twopi * real(iconnect) + gzz(iz), gvl(iv),           &
                  abs( conjg( ck(my)**iconnect * gff(-mxw,-my,iz,iv,im,1:ns) ) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( ofzv, "(99G17.7E3)" )                               &
                  + twopi * real(iconnect) + gzz(iz), gvl(iv),           &
                  abs( conjg( ck(my)**iconnect ) * gff(mxw,my,iz,iv,im,1:ns) )
              end do

          end if

        end do

      write( ofzv, * )
    end do

    end if

    close( ofzv )

  end if                           !%%% to restart

END SUBROUTINE ffinzv


SUBROUTINE analysis_zvm( inum )
!-------------------------------------------------------------------------------
!
!     Data analysis from gkv.****.zvm.%%%
!
!                                   by S.Maeyama  (Aug. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer, intent(in) :: inum

! --- local variables

  character*3 :: cnum
  integer :: ir, irz, irv, irm   ! rank, rankz, rankv, rankm
  character*4 :: crank
  integer :: loop, ios

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(-nz:nz-1,1:2*nv,0:nm,1:ns,1:15)         :: fzvm
  complex(kind=DP),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm,1:ns,1:15) :: gfzvm
  integer :: iz, iv, im, is, iw, giz, giv, gim

    write( cnum, fmt="(i3.3)" ) inum
    write( olog, * ) " # Analysis of zvm."

! --- file open ---
    do ir = 0, nprocz * nprocv * nprocm - 1
      write( crank, fmt="(i4.4)" ) ir
      open( unit=1000+ir, file="./bdata/gkv."//crank//".zvm."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".zvm."//cnum, 1000+ir
    end do


! --- time-step loop ---
    loop = 0
    do

  if ( loop >= 0 ) then            !%%% to restart

      do ir = 0, nprocz * nprocv * nprocm - 1

        read( unit=1000+ir, iostat=ios ) time, fzvm
        if ( ios /= 0 ) exit

        irz = mod( ir, nprocz )
        irv = mod( ir / nprocz, nprocv )
        irm = int( ir / ( nprocz * nprocv ) )
        do iw = 1, 15
          do is = 1, ns
            do im = 0, nm
              gim = (nm+1) * irm + im
              do iv = 1, 2*nv
                giv = 2*nv * irv + iv
                do iz = -nz, nz-1
                  giz = - global_nz + 2*nz * irz + iz + nz
                  gfzvm(giz,giv,gim,is,iw) = fzvm(iz,iv,im,is,iw)
                end do
              end do
            end do
          end do
        end do

      end do

  else                             !%%% to restart

      do ir = 0, nprocz * nprocv * nprocm - 1

        read( unit=1000+ir, iostat=ios )
        if ( ios /= 0 ) exit

      end do

  end if                           !%%% to restart

      if ( ios < 0 ) then
        write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                         ios, 1000+ir
        write( olog, * ) ""
        exit
      else if ( ios > 0 ) then
        write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                         ios, 1000+ir
        write( olog, * ) ""
        exit
      end if

      call fzvminzv( inum, loop, time, gfzvm )
      call fzvminvm( inum, loop, time, gfzvm )
      write( olog, * ) " # Data output at loop, time = ", loop, time

      loop = loop + 1

    end do


! --- file close ---
    do ir = 0, nprocz * nprocv * nprocm - 1
      close( unit=1000+ir )
    end do


END SUBROUTINE analysis_zvm


SUBROUTINE fzvminzv( inum, loop, time, gfzvm )
!-------------------------------------------------------------------------------
!
!     Write ff in zz, vl space
!
!                                   by S.Maeyama  (Aug. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                       :: inum
  integer,          intent(in)                       :: loop
  real(kind=DP),    intent(in)                       :: time
  complex(kind=DP), intent(in),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm,1:ns,1:15) :: gfzvm

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: iz, iv, im, iw

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pzzv, * ) "# set output './plt/fzvminzv_t"//cloop//"."//cnum//".jpg'"
    write( pzzv, * ) "splot './plt/fzvminzv_t"//cloop//"."//cnum//"' u 1:2:3 w pm3d"
    write( pzzv, * ) "pause 0.2"

    open( ozzv, file="./plt/fzvminzv_t"//cloop//"."//cnum )
    write( ozzv, * ) "#  time = ", time 
    write( ozzv, "(99a17)" ) "#              zz","vl","|ff|"
    im = 4
    iw = 1
      do iv = 1, 2*global_nv
        do iz = -global_nz, global_nz-1
          write( ozzv, "(99G17.7E3)" ) gzz(iz), gvl(iv),  &
                                       abs( gfzvm(iz,iv,im,1:ns,iw) )
        end do
        write( ozzv, * )
      end do
    close( ozzv )


END SUBROUTINE fzvminzv


SUBROUTINE fzvminvm( inum, loop, time, gfzvm )
!-------------------------------------------------------------------------------
!
!     Write ff in zz, vl space
!
!                                   by S.Maeyama  (Aug. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                       :: inum
  integer,          intent(in)                       :: loop
  real(kind=DP),    intent(in)                       :: time
  complex(kind=DP), intent(in),  &
    dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm,1:ns,1:15) :: gfzvm

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: iz, iv, im, iw

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pzvm, * ) "# set output './plt/fzvminvm_t"//cloop//"."//cnum//".jpg'"
    write( pzvm, * ) "splot './plt/fzvminvm_t"//cloop//"."//cnum//"' u 1:2:3 w pm3d"
    write( pzvm, * ) "pause 0.2"

    open( ozvm, file="./plt/fzvminvm_t"//cloop//"."//cnum )
    write( ozvm, * ) "#  time = ", time 
    write( ozvm, "(99a17)" ) "#              vl","mu","|ff|"
    iz = 0
    iw = 1
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          write( ozvm, "(99G17.7E3)" ) gvl(iv), gmu(im),  &
                                       abs( gfzvm(iz,iv,im,1:ns,iw) )
        end do
        write( ozvm, * )
      end do
    close( ozvm )


END SUBROUTINE fzvminvm


SUBROUTINE analysis_mom( inum )
!-------------------------------------------------------------------------------
!
!     Data analysis from gkv.****.phi.%%%
!     Data analysis from gkv.****.Al.%%%
!     Data analysis from gkv.****.mom.%%%
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer, intent(in) :: inum

! --- local variables

  character*3 :: cnum
  integer :: ir, irz             ! rank, rankz
  character*4 :: crank
  integer :: loop, ios

  real(kind=DP) :: time
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-nz:nz-1)               :: phi, Al
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-nz:nz-1,1:ns)          :: dens, upara,  &
                                                     ppara, pperp, &
                                                     qlpara, qlperp
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi, gAl
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1,1:ns) :: gdens, gupara,  &
                                                          gppara, gpperp, &
                                                          gqlpara, gqlperp
  integer :: mx, my, iz, giz, is

    write( cnum, fmt="(i3.3)" ) inum
    write( olog, * ) " # Analysis mom."

! --- file open ---
    do ir = 0, nprocz - 1
      write( crank, fmt="(i4.4)" ) ir
      open( unit=1000+ir, file="./bdata/gkv."//crank//".phi."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".phi."//cnum, 1000+ir
      open( unit=2000+ir, file="./bdata/gkv."//crank//".Al."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".Al."//cnum, 2000+ir
      open( unit=3000+ir, file="./bdata/gkv."//crank//".mom."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".mom."//cnum, 3000+ir
    end do


! --- time-step loop ---
    loop = 0
    do

  if ( loop >= 0 ) then            !%%% to restart

      do ir = 0, nprocz - 1

        read( unit=1000+ir, iostat=ios ) time, phi
        read( unit=2000+ir, iostat=ios ) time, Al
        read( unit=3000+ir, iostat=ios ) time, dens, upara, ppara, pperp, qlpara, qlperp
        if ( ios /= 0 ) exit

        irz = mod( ir, nprocz )
        do iz = -nz, nz-1
          giz = - global_nz + 2*nz * irz + iz + nz
          do my = -ny, ny
            do mx = 0, nx
              gphi(mx,my,giz) = phi(mx,my,iz)
               gAl(mx,my,giz) =  Al(mx,my,iz)
            end do
          end do
        end do
        do is = 1, ns
          do iz = -nz, nz-1
            giz = - global_nz + 2*nz * irz + iz + nz
            do my = -ny, ny
              do mx = 0, nx
                  gdens(mx,my,giz,is) =   dens(mx,my,iz,is)
                 gupara(mx,my,giz,is) =  upara(mx,my,iz,is)
                 gppara(mx,my,giz,is) =  ppara(mx,my,iz,is)
                 gpperp(mx,my,giz,is) =  pperp(mx,my,iz,is)
                gqlpara(mx,my,giz,is) = qlpara(mx,my,iz,is)
                gqlperp(mx,my,giz,is) = qlperp(mx,my,iz,is)
              end do
            end do
          end do
        end do

      end do

  else                             !%%% to restart

      do ir = 0, nprocz - 1

        read( unit=1000+ir, iostat=ios )
        if ( ios /= 0 ) exit
        read( unit=2000+ir, iostat=ios )
        if ( ios /= 0 ) exit
        read( unit=3000+ir, iostat=ios )
        if ( ios /= 0 ) exit

      end do

  end if                           !%%% to restart

      if ( ios < 0 ) then
        write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                         ios, 1000+ir, 2000+ir, 3000+ir
        write( olog, * ) ""
        exit
      else if ( ios > 0 ) then
        write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                         ios, 1000+ir, 2000+ir, 3000+ir
        write( olog, * ) ""
        exit
      end if

      call phiinkxky( inum, loop, time, gphi )
      call phiinxy( inum, loop, time, gphi )
      call phiinzz( inum, loop, time, gphi )
      call Alinkxky( inum, loop, time, gAl )
      call Alinxy( inum, loop, time, gAl )
      call Alinzz( inum, loop, time, gAl )

      call fluxinkxky( inum, loop, time, gphi, gAl,  &
                       gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
      call fluxintx( inum, loop, time, gphi, gAl,  &
                     gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
      write( olog, * ) " # Data output at loop, time = ", loop, time

      loop = loop + 1

    end do


! --- file close ---
    do ir = 0, nprocz - 1
      close( unit=1000+ir )
    end do


END SUBROUTINE analysis_mom


SUBROUTINE phiinkxky( inum, loop, time, gphi )
!-------------------------------------------------------------------------------
!
!     Write phi in kx, ky space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( ppkk, * ) "# set output './plt/phiinkxky_t"//cloop//"."//cnum//".jpg'"
    write( ppkk, * ) "splot './plt/phiinkxky_t"//cloop//"."//cnum//"' u 1:2:3"
    write( ppkk, * ) "pause 0.2"

    open( opkk, file="./plt/phiinkxky_t"//cloop//"."//cnum )
    write( opkk, * ) "#  time = ", time 
    write( opkk, "(99a17)" ) "#              kx","ky","|phi|"
    iz = 0
    do my = -ny ,ny
      do mx = 0, nx
        write( opkk, "(99G17.7E3)" ) kx(mx), ky(my), abs( gphi(mx,my,iz) )
      end do
      write( opkk, * )
    end do
    close( opkk )


END SUBROUTINE phiinkxky


SUBROUTINE phiinxy( inum, loop, time, gphi )
!-------------------------------------------------------------------------------
!
!     Write phi in x, y space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header
  use analysis_fft, only : fft_backward

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi

! --- local variables

  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: phixy
  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz


    do my = 0, 2*nyw-1
      do mx = 0, 2*nxw-1
        phixy(mx,my) = 0._DP
      end do
    end do

    iz = 0
    call fft_backward( gphi(:,:,iz), phixy )

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( ppxy, * ) "# set output './plt/phiinxy_t"//cloop//"."//cnum//".jpg'"
    write( ppxy, * ) "splot './plt/phiinxy_t"//cloop//"."//cnum//"' u 1:2:3"
    write( ppxy, * ) "pause 0.2"

    open( opxy, file="./plt/phiinxy_t"//cloop//"."//cnum )
    write( opxy, * ) "#  time = ", time 
    write( opxy, "(99a17)" ) "#              x","y","phixy"
    do my = 0, 2*nyw-1
      do mx = 0, 2*nxw-1
        write( opxy, "(99G17.7E3)" ) xx(mx), yy(my), phixy(mx,my)
      end do
      write( opxy, * )
    end do
    close( opxy )


END SUBROUTINE phiinxy


SUBROUTINE phiinzz( inum, loop, time, gphi )
!-------------------------------------------------------------------------------
!
!     Write phi in zz space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: mx, my, iz

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( ppzz, * ) "# set output './plt/phiinzz_t"//cloop//"."//cnum//".jpg'"
    write( ppzz, * ) "plot './plt/phiinzz_t"//cloop//"."//cnum//"' u 3:4 w l, '' u 3:5 w l"
    write( ppzz, * ) "pause 0.2"

  if ( loop >= 0 ) then            !%%% to restart

    open( opzz, file="./plt/phiinzz_t"//cloop//"."//cnum )
    write( opzz, * ) "#  time = ", time 
    write( opzz, "(99a17)" ) "#              kx","ky","zz","Re[phi]","Im[phi]"
    mx = 0
    my = 4

    if ( dj(my) == 0 ) then

        do iz = -global_nz, global_nz-1
          write( opzz, "(99G17.7E3)" ) kx(mx), ky(my), gzz(iz),  &
                                     gphi(mx,my,iz)
        end do

    else

      connect_min = int( ( nx + mx ) / abs( dj(my) ) )
      if ( connect_min .ne. 0 ) then
        do iconnect = connect_min, 1, -1
          mxw = mx+iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( opzz, "(99G17.7E3)" )                                &
                  kx(-mxw), ky(-my), - twopi * real(iconnect) + gzz(iz),  &
                  ck(my)**iconnect * conjg( gphi(-mxw,-my,iz) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( opzz, "(99G17.7E3)" )                              &
                  kx(mxw), ky(my), - twopi * real(iconnect) + gzz(iz),  &
                  ck(my)**iconnect * gphi(mxw,my,iz)
              end do

          end if

        end do
      end if

      connect_max = int( ( nx - mx ) / abs( dj(my) ) )
        do iconnect = 0, connect_max
          mxw = mx-iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( opzz, "(99G17.7E3)" )                                &
                  kx(-mxw), ky(-my), + twopi * real(iconnect) + gzz(iz),  &
                  conjg( ck(my)**iconnect * gphi(-mxw,-my,iz) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( opzz, "(99G17.7E3)" )                              &
                  kx(mxw), ky(my), + twopi * real(iconnect) + gzz(iz),  &
                  conjg( ck(my)**iconnect ) * gphi(mxw,my,iz)
              end do

          end if

        end do

    end if

    close( opzz )

  end if                           !%%% to restart


END SUBROUTINE phiinzz


SUBROUTINE Alinkxky( inum, loop, time, gAl )
!-------------------------------------------------------------------------------
!
!     Write Al in kx, ky space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gAl

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pakk, * ) "# set output './plt/Alinkxky_t"//cloop//"."//cnum//".jpg'"
    write( pakk, * ) "splot './plt/Alinkxky_t"//cloop//"."//cnum//"' u 1:2:3"
    write( pakk, * ) "pause 0.2"

    open( oakk, file="./plt/Alinkxky_t"//cloop//"."//cnum )
    write( oakk, * ) "#  time = ", time 
    write( oakk, "(99a17)" ) "#              kx","ky","|Al|"
    iz = 0
    do my = -ny ,ny
      do mx = 0, nx
        write( oakk, "(99G17.7E3)" ) kx(mx), ky(my), abs( gAl(mx,my,iz) )
      end do
      write( oakk, * )
    end do
    close( oakk )


END SUBROUTINE Alinkxky


SUBROUTINE Alinxy( inum, loop, time, gAl )
!-------------------------------------------------------------------------------
!
!     Write Al in x, y space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header
  use analysis_fft, only : fft_backward

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gAl

! --- local variables

  real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: Alxy
  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz


    do my = 0, 2*nyw-1
      do mx = 0, 2*nxw-1
        Alxy(mx,my) = 0._DP
      end do
    end do

    iz = 0
    call fft_backward( gAl(:,:,iz), Alxy )

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( paxy, * ) "# set output './plt/Alinxy_t"//cloop//"."//cnum//".jpg'"
    write( paxy, * ) "splot './plt/Alinxy_t"//cloop//"."//cnum//"' u 1:2:3"
    write( paxy, * ) "pause 0.2"

    open( oaxy, file="./plt/Alinxy_t"//cloop//"."//cnum )
    write( oaxy, * ) "#  time = ", time 
    write( oaxy, "(99a17)" ) "#              x","y","Alxy"
    do my = 0, 2*nyw-1
      do mx = 0, 2*nxw-1
        write( oaxy, "(99G17.7E3)" ) xx(mx), yy(my), Alxy(mx,my)
      end do
      write( oaxy, * )
    end do
    close( oaxy )


END SUBROUTINE Alinxy


SUBROUTINE Alinzz( inum, loop, time, gAl )
!-------------------------------------------------------------------------------
!
!     Write Al in zz space
!
!                                   by S.Maeyama  (Sep. 2010)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)                    :: inum
  integer,          intent(in)                    :: loop
  real(kind=DP),    intent(in)                    :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gAl

! --- local variables

  character*3 :: cnum
  character*8 :: cloop
  integer :: iconnect, connect_min, connect_max, mxw
  integer :: mx, my, iz

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pazz, * ) "# set output './plt/Alinzz_t"//cloop//"."//cnum//".jpg'"
    write( pazz, * ) "plot './plt/Alinzz_t"//cloop//"."//cnum//"' u 3:4 w l, '' u 3:5 w l"
    write( pazz, * ) "pause 0.2"

  if ( loop >= 0 ) then            !%%% to restart

    open( oazz, file="./plt/Alinzz_t"//cloop//"."//cnum )
    write( oazz, * ) "#  time = ", time 
    write( oazz, "(99a17)" ) "#              kx","ky","zz","Re[Al]","Im[Al]"
    mx = 0
    my = 4

    if ( dj(my) == 0 ) then

        do iz = -global_nz, global_nz-1
          write( oazz, "(99G17.7E3)" ) kx(mx), ky(my), gzz(iz),  &
                                     gAl(mx,my,iz)
        end do

    else

      connect_min = int( ( nx + mx ) / abs( dj(my) ) )
      if ( connect_min .ne. 0 ) then
        do iconnect = connect_min, 1, -1
          mxw = mx+iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( oazz, "(99G17.7E3)" )                                &
                  kx(-mxw), ky(-my), - twopi * real(iconnect) + gzz(iz),  &
                  ck(my)**iconnect * conjg( gAl(-mxw,-my,iz) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( oazz, "(99G17.7E3)" )                              &
                  kx(mxw), ky(my), - twopi * real(iconnect) + gzz(iz),  &
                  ck(my)**iconnect * gAl(mxw,my,iz)
              end do

          end if

        end do
      end if

      connect_max = int( ( nx - mx ) / abs( dj(my) ) )
        do iconnect = 0, connect_max
          mxw = mx-iconnect*dj(my)

          if ( mxw < 0 ) then

              do iz = -global_nz, global_nz-1
                write( oazz, "(99G17.7E3)" )                                &
                  kx(-mxw), ky(-my), + twopi * real(iconnect) + gzz(iz),  &
                  conjg( ck(my)**iconnect * gAl(-mxw,-my,iz) )
              end do

          else

              do iz = -global_nz, global_nz-1
                write( oazz, "(99G17.7E3)" )                              &
                  kx(mxw), ky(my), + twopi * real(iconnect) + gzz(iz),  &
                  conjg( ck(my)**iconnect ) * gAl(mxw,my,iz)
              end do

          end if

        end do

    end if

    close( oazz )

  end if                           !%%% to restart


END SUBROUTINE Alinzz


SUBROUTINE fluxinkxky( inum, loop, time, gphi, gAl,  &
                       gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!-------------------------------------------------------------------------------
!
!     Write particle and heat fluxes in kx, ky
!
!                                   by S.Maeyama  (Aug. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- argument

  integer,          intent(in)  :: inum
  integer,          intent(in)  :: loop
  real(kind=DP),    intent(in)  :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1)      :: gphi, gAl
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1,1:ns) :: gdens, gupara,  &
                                                          gppara, gpperp, &
                                                          gqlpara, gqlperp

! --- local variables

  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,1:ns) :: p_flux_es, p_flux_em, h_flux_es, h_flux_em
  real(kind=DP),  &
    dimension(1:ns) :: tp_flux_es, tp_flux_em, th_flux_es, th_flux_em
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: wc3
  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz, is


    tp_flux_es(:) = 0._DP
    tp_flux_em(:) = 0._DP
    th_flux_es(:) = 0._DP
    th_flux_em(:) = 0._DP

!%%% Electrostatic part of particle fluxes %%%
    do is = 1, ns
      do iz = -global_nz, global_nz-1
        do my = -ny, ny
          do mx = 0, nx
            wc3(mx,my,iz) = conjg( - ui * ky(my) * gphi(mx,my,iz) ) &
                                                 * gdens(mx,my,iz,is)
          end do
        end do
      end do
      call thet_ave_z ( wc3, p_flux_es(:,:,is) )
      do my = -ny, ny
        do mx = 1, nx
          tp_flux_es(is) = tp_flux_es(is) + T_s(is) * ( 1._DP + eta_s(is) )  &
                                * 2._DP * real( p_flux_es(mx,my,is), kind=DP )
        end do
      end do
      mx   = 0
        do my = 1, ny
          tp_flux_es(is) = tp_flux_es(is) + T_s(is) * ( 1._DP + eta_s(is) )  &
                                * 2._DP * real( p_flux_es(mx,my,is), kind=DP )
        end do
    end do


!%%% Magnetic part of particle fluxes %%%
    do is = 1, ns
      do iz = -global_nz, global_nz-1
        do my = -ny, ny
          do mx = 0, nx
            wc3(mx,my,iz) = conjg(   ui * ky(my) * gAl(mx,my,iz) ) &
                                               * gupara(mx,my,iz,is)
          end do
        end do
      end do
      call thet_ave_z ( wc3, p_flux_em(:,:,is) )
      do my = -ny, ny
        do mx = 1, nx
          tp_flux_em(is) = tp_flux_em(is) + T_s(is) * ( 1._DP + eta_s(is) )  &
                                * 2._DP * real( p_flux_em(mx,my,is), kind=DP )
        end do
      end do
      mx   = 0
        do my = 1, ny
          tp_flux_em(is) = tp_flux_em(is) + T_s(is) * ( 1._DP + eta_s(is) )  &
                                * 2._DP * real( p_flux_em(mx,my,is), kind=DP )
        end do
    end do


!%%% Electrostatic part of heat fluxes %%%
    do is = 1, ns
      do iz = -global_nz, global_nz-1
        do my = -ny, ny
          do mx = 0, nx
            wc3(mx,my,iz) = conjg( - ui * ky(my) * gphi(mx,my,iz) ) &
                      * ( gppara(mx,my,iz,is) + gpperp(mx,my,iz,is) &
                            - 2.5_DP * T_s(is) * gdens(mx,my,iz,is) )
          end do
        end do
      end do
      call thet_ave_z ( wc3, h_flux_es(:,:,is) )
      do my = -ny, ny
        do mx = 1, nx
          th_flux_es(is) = th_flux_es(is) + eta_s(is)  &
                                * 2._DP * real( h_flux_es(mx,my,is), kind=DP )
        end do
      end do
      mx   = 0
        do my = 1, ny
          th_flux_es(is) = th_flux_es(is) + eta_s(is)  &
                                * 2._DP * real( h_flux_es(mx,my,is), kind=DP )
        end do
    end do


!%%% Magnetic part of heat fluxes %%%
    do is = 1, ns
      do iz = -global_nz, global_nz-1
        do my = -ny, ny
          do mx = 0, nx
            wc3(mx,my,iz) = conjg(   ui * ky(my) * gAl(mx,my,iz) ) &
                   * ( gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is) &
                          - 2.5_DP * T_s(is) * gupara(mx,my,iz,is) )
          end do
        end do
      end do
      call thet_ave_z ( wc3, h_flux_em(:,:,is) )
      do my = -ny, ny
        do mx = 1, nx
          th_flux_em(is) = th_flux_em(is) + eta_s(is)  &
                                * 2._DP * real( h_flux_em(mx,my,is), kind=DP )
        end do
      end do
      mx   = 0
        do my = 1, ny
          th_flux_em(is) = th_flux_em(is) + eta_s(is)  &
                                * 2._DP * real( h_flux_em(mx,my,is), kind=DP )
        end do
    end do

    write( oxtk, "(99G17.7E3)" ) time, tp_flux_es(1:ns), tp_flux_em(1:ns),  &
                                       th_flux_es(1:ns), th_flux_em(1:ns)


!%%% Spectra of fluxes %%%
    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

    write( pxkk, * ) "# set output './plt/fluxinkxky_t"//cloop//"."//cnum//".jpg'"
    write( pxkk, * ) "splot './plt/fluxinkxky_t"//cloop//"."//cnum//"' u 1:2:3"
    write( pxkk, * ) "pause 0.2"

    open( oxkk, file="./plt/fluxinkxky_t"//cloop//"."//cnum )
    write( oxkk, * ) "#  time = ", time 
    write( oxkk, "(99a17)" ) "#              kx","ky",            &
                             "|p_flux_es(ns)|","|p_flux_em(ns)|", &
                             "|h_flux_es(ns)|","|h_flux_em(ns)|"
    do my = -ny ,ny
      do mx = 0, nx
        write( oxkk, "(99G17.7E3)" ) kx(mx), ky(my),                         &
                abs( p_flux_es(mx,my,1:ns) ), abs( p_flux_em(mx,my,1:ns) ),  &
                abs( h_flux_es(mx,my,1:ns) ), abs( h_flux_em(mx,my,1:ns) )
      end do
      write( oxkk, * )
    end do
    close( oxkk )


END SUBROUTINE fluxinkxky


SUBROUTINE fluxintx( inum, loop, time, gphi, gAl,  &
                     gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!-------------------------------------------------------------------------------
!
!     Write particle and heat fluxes in x, y
!
!                                   by S.Maeyama  (Aug. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header
  use analysis_fft, only : fft_backward

  implicit none

! --- argument

  integer,          intent(in)  :: inum
  integer,          intent(in)  :: loop
  real(kind=DP),    intent(in)  :: time
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1)      :: gphi, gAl
  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1,1:ns) :: gdens, gupara,  &
                                                          gppara, gpperp, &
                                                          gqlpara, gqlperp

! --- local variables

  real(kind=DP),  &
    dimension(0:2*nxw-1,0:2*nyw-1) :: exbxy, mfltxy
  real(kind=DP),  &
    dimension(0:2*nxw-1,0:2*nyw-1,1:ns) :: densxy, uparaxy,  &
                                           pparaxy, pperpxy, &
                                           qlparaxy, qlperpxy
  real(kind=DP),  &
    dimension(0:2*nxw-1,1:ns) :: p_flux_es, p_flux_em,  &
                                 h_flux_es, h_flux_em
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny) :: wc2
  character*3 :: cnum
  character*8 :: cloop
  integer :: mx, my, iz, is


    iz = 0
!%%% ExB flow %%%
      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = - ui * ky(my) * gphi(mx,my,iz)
        end do
      end do
      call fft_backward( wc2, exbxy )


!%%% Magnetic flutter %%%
      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = ui * ky(my) * gAl(mx,my,iz)
        end do
      end do
      call fft_backward( wc2, mfltxy )


!%%% Fluctuations of density, flow, pressure, heat flux %%%
    do is = 1, ns

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gdens(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, densxy(:,:,is) )

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gupara(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, uparaxy(:,:,is) )

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gppara(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, pparaxy(:,:,is) )

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gpperp(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, pperpxy(:,:,is) )

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gqlpara(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, qlparaxy(:,:,is) )

      do my = -ny, ny
        do mx = 0, nx
          wc2(mx,my) = gqlperp(mx,my,iz,is)
        end do
      end do
      call fft_backward( wc2, qlperpxy(:,:,is) )

    end do


!%%% Particle and heat fluxes %%%
    p_flux_es(:,:) = 0._DP
    p_flux_em(:,:) = 0._DP
    h_flux_es(:,:) = 0._DP
    h_flux_em(:,:) = 0._DP
    do is = 1, ns
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          p_flux_es(mx,is) = p_flux_es(mx,is) + exbxy(mx,my) * densxy(mx,my,is)
          p_flux_em(mx,is) = p_flux_em(mx,is) + mfltxy(mx,my) * uparaxy(mx,my,is)
          h_flux_es(mx,is) = h_flux_es(mx,is) + exbxy(mx,my) &
                                    * ( pparaxy(mx,my,is) + pperpxy(mx,my,is) &
                                        - 2.5_DP * T_s(is) * densxy(mx,my,is) )
          h_flux_em(mx,is) = h_flux_em(mx,is) + mfltxy(mx,my) &
                                    * ( qlparaxy(mx,my,is) + qlperpxy(mx,my,is) &
                                         - 2.5_DP * T_s(is) * uparaxy(mx,my,is) )
        end do
      end do
    end do
    p_flux_es(:,:) = p_flux_es(:,:) / real( 2*nyw, kind=DP )
    p_flux_em(:,:) = p_flux_em(:,:) / real( 2*nyw, kind=DP )
    h_flux_es(:,:) = h_flux_es(:,:) / real( 2*nyw, kind=DP )
    h_flux_em(:,:) = h_flux_em(:,:) / real( 2*nyw, kind=DP )

    do mx = 0, 2*nxw-1
      write( oxtx, "(99G17.7E3)" ) time, xx(mx),              &
                     p_flux_es(mx,1:ns), p_flux_em(mx,1:ns),  &
                     h_flux_es(mx,1:ns), h_flux_em(mx,1:ns)
    end do
    write( oxtx, * )


END SUBROUTINE fluxintx


SUBROUTINE thet_ave_z ( wn, wa )
!-------------------------------------------------------------------------------
!
!     average of a complex variable wn in the theta space
!
!-------------------------------------------------------------------------------

  use analysis_header

  implicit none

! --- arguments

  complex(kind=DP), intent(in),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: wn

  complex(kind=DP), intent(out), &
    dimension(0:nx,-ny:ny)                        :: wa

! --- local variables

  real(kind=DP) :: fct
  integer ::  mx, my, iz


    wa   = ( 0._DP, 0._DP )

    do iz = -global_nz, global_nz-1
      fct = 1._DP / ( cfsrf * gbeq(iz) )
      do my = -ny, ny
        do mx = 0, nx
          wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
        end do
      end do
    end do


END SUBROUTINE thet_ave_z


END PROGRAM analysis
