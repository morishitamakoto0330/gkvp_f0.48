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

! --- unit numbers for I/O
  integer, parameter :: olog = 10, &
                        iset = 12, &
                        ofzv = 20, pfzv = 21, &
                        opkk = 30, ppkk = 31, &
                        opxy = 32, ppxy = 33, &
                        opzz = 34, ppzz = 35

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

    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: phikk
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: phixy


      call dfftw_plan_dft_c2r_2d( plan, 2*nxw, 2*nyw, &
                                  phikk, phixy,       &
                                  FFTW_ESTIMATE )


  END SUBROUTINE fft_pre


!--------------------------------------
  SUBROUTINE fft_backward ( gphi, phixy )
!--------------------------------------
!  Execution of FFT

    complex(kind=DP), intent(in), &
      dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi
    real(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:2*nyw-1)                :: phixy

    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: phikk
    integer :: mx, my, iz

    iz = 0
      do my = 0, 2*nyw-1
        do mx = 0, nxw
          phikk(mx,my) = ( 0._DP, 0._DP )
        end do
      end do
      do my = 0, ny
        do mx = 0, nx
          phikk(mx,my) = gphi(mx,my,iz)
        end do
      end do
      do my = -ny, -1
        do mx = 0, nx
          phikk(mx,2*nyw+my) = gphi(mx,my,iz)
        end do
      end do
      mx = 0
        do my = 1, ny
          phikk(mx,2*nyw-my) = conjg( phikk(mx,my) )
        end do

    call dfftw_execute_dft( plan, phikk, phixy )


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
      call analysis_phi( inum )

      if ( inum == enum ) then
        call reset_parameters
      end if

      write( olog, * ) " ######################"
      write( olog, * ) ""

    end do

    call file_close


END PROGRAM analysis


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
    close( ppkk )
    close( ppxy )
    close( ppzz )


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

    write( pfzv, * ) "splot './plt/ffinzv_t"//cloop//"."//cnum//"' u 1:2:3 w pm3d"
    write( pfzv, * ) "pause 0.2"

  if ( loop >= 0 ) then            !%%% to restart

    open( ofzv, file="./plt/ffinzv_t"//cloop//"."//cnum )
    write( ofzv, * ) "#  time = ", time 
    write( ofzv, "(99a17)" ) "#              zz","vl","|ff|"
    mx = 0
    my = 1
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


SUBROUTINE analysis_phi( inum )
!-------------------------------------------------------------------------------
!
!     Data analysis from gkv.****.phi.%%%
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
    dimension(0:nx,-ny:ny,-nz:nz-1)               :: phi
  complex(kind=DP),  &
    dimension(0:nx,-ny:ny,-global_nz:global_nz-1) :: gphi
  integer :: mx, my, iz, giz

    write( cnum, fmt="(i3.3)" ) inum
    write( olog, * ) " # Analysis phi."

! --- file open ---
    do ir = 0, nprocz - 1
      write( crank, fmt="(i4.4)" ) ir
      open( unit=1000+ir, file="./bdata/gkv."//crank//".phi."//cnum,  &
            status="old", action="read", form="unformatted" )
      write( olog, * ) " # Opened file and the unit-ID are ",  &
                       "./bdata/gkv."//crank//".phi."//cnum, 1000+ir
    end do


! --- time-step loop ---
    loop = 0
    do

  if ( loop >= 0 ) then            !%%% to restart

      do ir = 0, nprocz - 1

        read( unit=1000+ir, iostat=ios ) time, phi
        if ( ios /= 0 ) exit

        irz = mod( ir, nprocz )
        do iz = -nz, nz-1
          giz = - global_nz + 2*nz * irz + iz + nz
          do my = -ny, ny
            do mx = 0, nx
              gphi(mx,my,giz) = phi(mx,my,iz)
            end do
          end do
        end do

      end do

  else                             !%%% to restart

      do ir = 0, nprocz - 1

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

      call phiinkxky( inum, loop, time, gphi )
      call phiinxy( inum, loop, time, gphi )
      call phiinzz( inum, loop, time, gphi )
      write( olog, * ) " # Data output at loop, time = ", loop, time

      loop = loop + 1

    end do


! --- file close ---
    do ir = 0, nprocz - 1
      close( unit=1000+ir )
    end do


END SUBROUTINE analysis_phi


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

    call fft_backward( gphi, phixy )

    write( cnum, fmt="(i3.3)" ) inum
    write( cloop, fmt="(i8.8)" ) loop

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

    write( ppzz, * ) "plot './plt/phiinzz_t"//cloop//"."//cnum//"' u 3:4 w l"
    write( ppzz, * ) "pause 0.2"

  if ( loop >= 0 ) then            !%%% to restart

    open( opzz, file="./plt/phiinzz_t"//cloop//"."//cnum )
    write( opzz, * ) "#  time = ", time 
    write( opzz, "(99a17)" ) "#              kx","ky","zz","|phi|"
    mx = 0
    my = 1

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


