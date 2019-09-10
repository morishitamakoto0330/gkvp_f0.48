MODULE analysis_header
!-------------------------------------------------------------------------------
!
!     Header for data analysis
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)


!%%% Set parameters %%%
  integer, parameter :: snum = 1      ! begining of simulation runs
  integer, parameter :: enum = 1      ! end of simulation runs

  integer, parameter :: nxw = 288, nyw = 48
  integer, parameter :: nx = 192, global_ny = 31  ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 32, global_nv = 32, global_nm = 15
  integer, parameter :: nzb = 2
  integer, parameter :: nprocw = 8, nprocz = 8, nprocv = 8, nprocm = 2, nprocs = 2
  real(kind=DP), parameter ::  vmax = 4._DP
!%%%%%%%%%%%%%%%%%%%%%%


  integer, parameter :: nxw_size = (2*nxw-1)/nprocw     ! local allocation size (0:nxw_size)
  integer, parameter :: ny       = global_ny / nprocw   ! local allocation size (0:ny)
  integer, parameter :: nz = global_nz / nprocz,          &
                        nv = global_nv / nprocv,          &
                        nm = (global_nm + 1) / nprocm - 1,&
                        ns = nprocs
  integer, parameter :: nproc = nprocw * nprocz * nprocv * nprocm * nprocs

  real(kind=DP), dimension(0:2*nxw-1)              :: xx
  real(kind=DP), dimension(0:2*nyw-1)              :: yy
  real(kind=DP), dimension(-nx:nx)                 :: kx
  real(kind=DP), dimension(0:global_ny)            :: gky
  real(kind=DP), dimension(-global_nz:global_nz-1) :: gzz, gomg
  real(kind=DP), dimension(1:2*global_nv)          :: gvl
  real(kind=DP), dimension(0:global_nm)            :: gmu
  real(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) :: gfmx

  complex(kind=DP), dimension(0:global_ny)         :: ck
  integer, dimension(0:global_ny)                  :: dj

  real(kind=DP), dimension(0:ns-1) ::     eta,  &    ! T-gradient
                                           nu,  &    ! collision freq.   
                                         Anum,  &    ! mass number
                                         Znum,  &    ! charge number     
                                          fcs,  &    ! charge-density fraction 
                                          sgn,  &    ! signs of charge   
                                          tau,  &    ! T-ratio
                                      r_major        ! R0/Lns 
  real(kind=DP) :: dpara, dv, cfsrf, lambda_i, beta
  character(8)  :: equib_type
  real(kind=DP) :: r_minor, q_0, q_d, s_hat
  real(kind=DP) :: eps_r
  real(kind=DP) :: rdeps00, eps_hor, lprd, mprd, lmmq, malpha
  real(kind=DP) :: eps_mor, eps_por, lprdm1, lprdp1, lmmqm1, lmmqp1
  real(kind=DP) :: eps_rnew, rdeps1_0, rdeps1_10, rdeps2_10, rdeps3_10
  real(kind=DP) :: lx, ly, lz, kxmin, kymin, dz, mmax, dm
  integer       :: n_alp, n_tht

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, twopi = pi * 2._DP
  real(kind=DP),    parameter :: eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )

! --- unit numbers for I/O
  integer, parameter :: olog = 10, &
                        opkk = 40, ppkk = 41, &
                        opxy = 42, ppxy = 43, &
                        opzz = 44, ppzz = 45, &
                        oakk = 50, pakk = 51, &
                        oaxy = 52, paxy = 53, &
                        oazz = 54, pazz = 55


END MODULE analysis_header


MODULE analysis_fft
!-------------------------------------------------------------------------------
!
!    FFT module using fftw3
!
!                                   by S.Maeyama  (Dec. 2012)
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
  
    complex(kind=DP), dimension(0:2*nxw-1,0:nyw) :: wwkk
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
      dimension(-nx:nx,0:global_ny)              :: gww
    real(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:2*nyw-1)             :: wwxy
  
    complex(kind=DP), dimension(0:2*nxw-1,0:nyw) :: wwkk
    integer :: mx, my
  
  
      do my = 0, nyw
        do mx = 0, 2*nxw-1
          wwkk(mx,my) = ( 0._DP, 0._DP )
        end do
      end do
      do my = 0, global_ny
        do mx = 0, nx
          wwkk(mx,my) = gww(mx,my)
        end do
      end do
      do my = 0, global_ny
        do mx = -nx, -1
          wwkk(2*nxw+mx,ny) = gww(mx,my)
        end do
      end do
      my = 0
        do mx = 1, nx
          wwkk(2*nxw-mx,my) = conjg( wwkk(mx,my) )
        end do
  
    !  call dfftw_execute_dft( plan, wwkk, wwxy )
      call dfftw_execute_dft_c2r( plan, wwkk, wwxy )
  
  
  END SUBROUTINE fft_backward


END MODULE analysis_fft


PROGRAM analysis
!-------------------------------------------------------------------------------
!
!     Data analysis from binary output
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  use analysis_header, only : snum, enum, olog
  use analysis_fft, only : fft_pre

  implicit none

  integer :: inum
  character*3 :: cnum

    call file_open

    do inum = snum, enum

      if ( inum == snum ) then
        call init( inum )
        call fft_pre
      end if

      write( cnum, fmt="(i3.3)" ) inum
      write( olog, * ) " ##### inum = "//cnum//" #####"

      call analysis_mom( inum )

      write( olog, * ) " ######################"
      write( olog, * ) ""

    end do

    call file_close

  stop

CONTAINS


!--------------------------------------
  SUBROUTINE file_open
!--------------------------------------
!     Open files
  
    use analysis_header
  
    implicit none
  
      open( olog, file="./plt/log.dat" )
  
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
  
  
  END SUBROUTINE file_open
  
  
!--------------------------------------
  SUBROUTINE file_close
!--------------------------------------
!     Close files
  
    use analysis_header
  
    implicit none
  
      close( olog )
      close( ppkk )
      close( ppxy )
      close( ppzz )
      close( pakk )
      close( paxy )
      close( pazz )
  
  
  END SUBROUTINE file_close
  
  
!--------------------------------------
  SUBROUTINE init( inum )
!--------------------------------------
!     Initial setting
  
    use analysis_header
  
    implicit none
  
    integer, intent(in) :: inum
  
    character*3 :: cnum
    integer :: mx, my, iz, iv, im
  
    namelist /equib/ equib_type
    namelist /physp/ eta,  &    ! T-gradient
                      nu,  &    ! collision freq.
                    Anum,  &    ! mass number
                    Znum,  &    ! charge number 
                     fcs,  &    ! charge-density fraction 
                     sgn,  &    ! signs of charge 
                     tau,  &    ! T-ratio Ts/T0, T0=reference temp.
                lambda_i,  &    ! Debye/rho 
                    beta        ! Beta value
  
    namelist /nperi/ n_alp, n_tht
    namelist /confp/ r_minor, r_major, eps_r, eps_rnew,     &
                     q_0, q_d, s_hat,                       &
                     lprd, mprd, eps_hor, eps_mor, eps_por, &
                     rdeps00, rdeps1_0, rdeps1_10,          &
                     rdeps2_10, rdeps3_10, malpha
  
  
      write( cnum, fmt="(i3.3)" ) inum
  
      open( 5, file="./gkvp_f0.26_namelist."//cnum, status="old", action="read" )
  
      read(5,nml=equib)
  
      read(5,nml=physp)
  
      read(5,nml=nperi)
  
      if( trim(equib_type) == "analytic" ) then
  
        read(5,nml=confp)
  
        lprdm1   = lprd - 1.0_DP
        lprdp1   = lprd + 1.0_DP
  
        lmmq     = lprd   - mprd * q_0
        lmmqm1   = lprdm1 - mprd * q_0
        lmmqp1   = lprdp1 - mprd * q_0
  
      else
  
        write( olog, * ) " # Currently, this is not available"
        stop
  
      end if
  
      lx       = r_minor * q_d / ( q_0 * s_hat )
      ly       = r_minor * pi  / ( q_0 * real( n_alp, kind=DP ) )
      lz       = real( n_tht, kind=DP ) * pi        ! total z-length
  
      kxmin    = pi / lx
      kymin    = pi / ly
  
      dz       = lz / real( global_nz, kind=DP )
      dpara    = dz * q_0 * r_major(0)
  
!      vmax     = 5._DP
      dv       = 2._DP * vmax / real( 2 * nv * nprocv -1, kind=DP )
  
      mmax     = vmax
      dm       = mmax / real( nprocm * ( nm+1 ) - 1, kind=DP )
  
      do mx = 0, 2*nxw-1
        xx(mx) = -lx + lx/real(nxw)*mx
      end do
  
      do my = 0, 2*nyw-1
        yy(my) = -ly + ly/real(nyw)*my
      end do
  
      do mx = -nx, nx
        kx(mx) = kxmin * real( mx, kind=DP )
      end do
  
      do my = 0, global_ny
        gky(my) = kymin * real( my, kind=DP )
      end do
  
      do iz = -global_nz, global_nz-1
        gzz(iz) = dz * real( iz, kind=DP )
        gomg(iz)  = 1._DP                                              &
                  - eps_r * ( cos( gzz(iz) )                           &
                          + eps_hor * cos( lmmq   * gzz(iz) - malpha ) &
                          + eps_mor * cos( lmmqm1 * gzz(iz) - malpha ) &
                          + eps_por * cos( lmmqp1 * gzz(iz) - malpha ) )
      end do
      cfsrf = 0._DP
      do iz = -global_nz, global_nz-1
        cfsrf = cfsrf + 1._DP / gomg(iz)
      end do
  
      do iv = 1, 2*global_nv
        gvl(iv) = - vmax + dv * real( iv-1, kind=DP )
      end do
  
      do im = 0, global_nm
        gmu(im) = 0.5_DP * ( dm * real( im, kind=DP ) )**2
      end do
  
      do my = 0, global_ny
        ck(my)   = exp( ui * 2._DP * pi * q_0 &
                           * real( n_alp * n_tht * my, kind=DP ) )
        dj(my)   = - nint( real( 2 * n_alp * n_tht * my, kind=DP ) * q_d )
      end do
  
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          do iz = -global_nz, global_nz-1
            gfmx(iz,iv,im) = exp( - 0.5_DP * gvl(iv)**2 - gomg(iz) * gmu(im) ) &
                           / sqrt( twopi**3 )
          end do
        end do
      end do
  
      close(5)
  
  
  END SUBROUTINE init
  
  
!--------------------------------------
  SUBROUTINE analysis_mom( inum )
!--------------------------------------
!     Data analysis from *.phi.*, *.Al.*, *.mom.*
  
    use analysis_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir, irw, irz
    character*6 :: crank
    character*1 :: srank
    integer :: loop, ios
  
    real(kind=DP) :: time
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor
    integer :: mx, my, gmy, iz, giz
  
      write( cnum, fmt="(i3.3)" ) inum
      write( olog, * ) " # Analysis mom."
  
  ! --- file open ---
      do ir = 0, nproc - 1
        rankg = ir
  
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
  
        scolor = mod( rankg, nprocz * nprocw )
  
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        
        if ( ranks == 0 .AND. rankv == 0 . AND. rankm == 0 ) then
  
          write( crank, fmt="(i6.6)" ) ir
          write( srank, fmt="(i1.1)" ) ranks
          open( unit=100000+ir, file="./phi/gkvp_f0.26."//crank//"."//srank//".phi."//cnum,  &
                status="old", action="read", form="unformatted" )
          write( olog, * ) " # Opened file and the unit-ID are ",  &
                           "./phi/gkvp_f0.26."//crank//"."//srank//".phi."//cnum, 100000+ir
          open( unit=200000+ir, file="./phi/gkvp_f0.26."//crank//"."//srank//".Al."//cnum,  &
                status="old", action="read", form="unformatted" )
          write( olog, * ) " # Opened file and the unit-ID are ",  &
                           "./phi/gkvp_f0.26."//crank//"."//srank//".Al."//cnum, 200000+ir
  
        end if
  
      end do
  
  
  ! --- time-step loop ---
      loop = 0
      do
  
!    if ( loop >= 0 ) then            !%%% to restart
    if ( loop <= 10 ) then            !%%% to restart
  
        do ir = 0, nproc - 1
          rankg = ir
    
          ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
          rank = mod( rankg, nproc / nprocs )
    
          scolor = mod( rankg, nprocz * nprocw )
    
          rankw = mod( rank,                      nprocw )
          rankz = mod( rank / nprocw,             nprocz )
          rankv = mod( rank / ( nprocw * nprocz), nprocv )
          rankm = rank / ( nprocw * nprocz * nprocv )
          wcolor   = rank / nprocw
          zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
          vcolor   = rankz * nprocw  +  rankw
          
          if ( ranks == 0 .AND. rankv == 0 .AND. rankm ==0 ) then
  
            read( unit=100000+ir, iostat=ios ) time, phi
            read( unit=200000+ir, iostat=ios ) time, Al
            if ( ios /= 0 ) exit
  
            do iz = -nz, nz-1
              giz = - global_nz + 2*nz * rankz + iz + nz
              do my = 0, ny
                gmy = ( ny+1 ) * rankw + my
                if ( gmy <= global_ny ) then
                  do mx = -nx, nx
                    gphi(mx,gmy,giz) = phi(mx,my,iz)
                     gAl(mx,gmy,giz) =  Al(mx,my,iz)
                  end do
                end if
              end do
            end do
  
          end if
  
        end do
  
    else                             !%%% to restart
  
        do ir = 0, nproc - 1
          rankg = ir
    
          ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
          rank = mod( rankg, nproc / nprocs )
    
          scolor = mod( rankg, nprocz * nprocw )
    
          rankw = mod( rank,                      nprocw )
          rankz = mod( rank / nprocw,             nprocz )
          rankv = mod( rank / ( nprocw * nprocz), nprocv )
          rankm = rank / ( nprocw * nprocz * nprocv )
          wcolor   = rank / nprocw
          zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
          vcolor   = rankz * nprocw  +  rankw
          
          if ( ranks == 0 .AND. rankv == 0 .AND. rankm ==0 ) then
  
            read( unit=100000+ir, iostat=ios )
            if ( ios /= 0 ) exit
            read( unit=200000+ir, iostat=ios )
            if ( ios /= 0 ) exit
  
          end if
  
        end do
  
    end if                           !%%% to restart
  
        if ( ios < 0 ) then
          write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                           ios, 100000+ir, 200000+ir
          write( olog, * ) ""
          exit
        else if ( ios > 0 ) then
          write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                           ios, 100000+ir, 200000+ir
          write( olog, * ) ""
          exit
        end if
  
!    if ( loop >= 0 ) then            !%%% to restart
    if ( loop <= 10 ) then            !%%% to restart

        call phiinkxky( inum, loop, time, gphi )
        call phiinxy( inum, loop, time, gphi )
        call phiinzz( inum, loop, time, gphi )
        call Alinkxky( inum, loop, time, gAl )
        call Alinxy( inum, loop, time, gAl )
        call Alinzz( inum, loop, time, gAl )
  
        write( olog, * ) " # Data output at loop, time = ", loop, time

    end if
  
        loop = loop + 1
  
      end do
  
  
  ! --- file close ---
      do ir = 0, nproc - 1
        rankg = ir
   
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
   
        scolor = mod( rankg, nprocz * nprocw )
   
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        
        if ( ranks == 0 .AND. rankv == 0 .AND. rankm ==0 ) then
  
          close( unit=100000+ir )
          close( unit=200000+ir )
  
        end if
  
      end do
  
  
  END SUBROUTINE analysis_mom
  
  
!--------------------------------------
  SUBROUTINE phiinkxky( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in kx, ky space
  
    use analysis_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
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
      do my = 0 ,global_ny
        do mx = -nx, nx
          write( opkk, "(99G17.7E3)" ) kx(mx), gky(my), abs( gphi(mx,my,iz) )
        end do
        write( opkk, * )
      end do
      close( opkk )
  
  
  END SUBROUTINE phiinkxky
  
  
!--------------------------------------
  SUBROUTINE phiinxy( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in x, y space
  
    use analysis_header
    use analysis_fft, only : fft_backward
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
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
  
  
!--------------------------------------
  SUBROUTINE phiinzz( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in zz space
  
    use analysis_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
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
  
      open( opzz, file="./plt/phiinzz_t"//cloop//"."//cnum )
      write( opzz, * ) "#  time = ", time 
      write( opzz, "(99a17)" ) "#              kx","ky","zz","Re[phi]","Im[phi]"
      mx = 0
      my = 4
  
      if ( dj(my) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( opzz, "(99G17.7E3)" ) kx(mx), gky(my), gzz(iz),  &
                                       gphi(mx,my,iz)
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(my) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( opzz, "(99G17.7E3)" )                              &
                kx(mxw), gky(my), - twopi * real(iconnect) + gzz(iz),  &
                ck(my)**iconnect * gphi(mxw,my,iz)
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(my) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( opzz, "(99G17.7E3)" )                              &
                kx(mxw), gky(my), + twopi * real(iconnect) + gzz(iz),  &
                conjg( ck(my)**iconnect ) * gphi(mxw,my,iz)
            end do
          end do
  
      end if
  
      close( opzz )
  
  
  END SUBROUTINE phiinzz
  
  
!--------------------------------------
  SUBROUTINE Alinkxky( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in kx, ky space
  
    use analysis_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
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
      do my = 0 ,global_ny
        do mx = -nx, nx
          write( oakk, "(99G17.7E3)" ) kx(mx), gky(my), abs( gAl(mx,my,iz) )
        end do
        write( oakk, * )
      end do
      close( oakk )
  
  
  END SUBROUTINE Alinkxky
  
  
!--------------------------------------
  SUBROUTINE Alinxy( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in x, y space
  
    use analysis_header
    use analysis_fft, only : fft_backward
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
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
  
  
!--------------------------------------
  SUBROUTINE Alinzz( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in zz space
  
    use analysis_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
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
  
      open( oazz, file="./plt/Alinzz_t"//cloop//"."//cnum )
      write( oazz, * ) "#  time = ", time 
      write( oazz, "(99a17)" ) "#              kx","ky","zz","Re[Al]","Im[Al]"
      mx = 0
      my = 4
  
      if ( dj(my) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( oazz, "(99G17.7E3)" ) kx(mx), gky(my), gzz(iz),  &
                                       gAl(mx,my,iz)
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(my) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( oazz, "(99G17.7E3)" )                              &
                kx(mxw), gky(my), - twopi * real(iconnect) + gzz(iz),  &
                ck(my)**iconnect * gAl(mxw,my,iz)
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(my) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( oazz, "(99G17.7E3)" )                              &
                kx(mxw), gky(my), + twopi * real(iconnect) + gzz(iz),  &
                conjg( ck(my)**iconnect ) * gAl(mxw,my,iz)
            end do
          end do
  
      end if
  
      close( oazz )
  
  
  END SUBROUTINE Alinzz


END PROGRAM analysis
