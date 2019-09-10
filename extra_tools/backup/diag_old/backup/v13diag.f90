MODULE diag_header
!-------------------------------------------------------------------------------
!
!     Header for data diagnostics
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  implicit none

  public

  integer, parameter :: DP = selected_real_kind(14)


!%%% DIAG parameters %%%
  integer, parameter :: snum = 1      ! begining of simulation runs
  integer, parameter :: enum = 1      ! end of simulation runs
  integer, parameter :: loopskip = 10
  logical, parameter :: flag_cnt = .false.
  logical, parameter :: flag_ffinzv    = .false.
  logical, parameter :: flag_mom = .true.
  logical, parameter :: flag_mominkxky = .false.
  logical, parameter :: flag_mominxy   = .true.
  logical, parameter :: flag_mominz    = .false.
  logical, parameter :: flag_trn = .true.
  logical, parameter :: flag_trninkxky = .true.
!%%%%%%%%%%%%%%%%%%%%%%

!%%% GKV parameters %%%
  integer, parameter :: nxw = 128, nyw = 32
  integer, parameter :: nx = 80, global_ny = 19  ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 32, global_nv = 48, global_nm = 15
  integer, parameter :: nzb = 2
  integer, parameter :: nprocw = 4, nprocz = 8, nprocv = 12, nprocm = 2, nprocs = 2
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
  real(kind=DP), dimension(-global_nz:global_nz-1) :: gzz, gomg, grootg
  real(kind=DP), dimension(1:2*global_nv)          :: gvl
  real(kind=DP), dimension(0:global_nm)            :: gmu
  real(kind=DP), dimension(-global_nz:global_nz-1,1:2*global_nv,0:global_nm) :: gfmx

  complex(kind=DP), dimension(0:global_ny)         :: ck
  integer, dimension(0:global_ny)                  :: dj

  real(kind=DP), dimension(0:ns-1) ::   R0_Ln,  &    ! R0/Lns
                                        R0_Lt,  &    ! R0/Lts
                                           nu,  &    ! collision freq.   
                                         Anum,  &    ! mass number
                                         Znum,  &    ! charge number     
                                          fcs,  &    ! charge-density fraction 
                                          sgn,  &    ! signs of charge   
                                          tau        ! T-ratio
  real(kind=DP) :: dv, cfsrf, lambda_i, beta, q_0, q_bar, theta
  character(8)  :: equib_type
  real(kind=DP) :: r_minor, q_d, s_hat
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
                     offinzv = 20, pffinzv = 21, &
                  omominkxky = 40, pmominkxky = 41, &
                    omominxy = 42, &
                     omominz = 44, &
                  otrninkxky = 70


END MODULE diag_header


MODULE diag_fft
!-------------------------------------------------------------------------------
!
!    FFT module using fftw3
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  use diag_header

  implicit none

  include "fftw3.f"

  private

  integer(kind=DP), save :: plan_backward_xy, plan_backward_x

  public   fft_pre, fft_backward_xy, fft_backward_x


CONTAINS


!--------------------------------------
  SUBROUTINE fft_pre
!--------------------------------------
!  Initialization of FFT
  
    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: wwxy

    complex(kind=DP), dimension(0:2*nxw-1) :: wwkx, wwxx

    integer :: ierr, omp_get_max_threads

!%%% setting for fftw3 %%%
!      call dfftw_plan_dft_c2r_2d( plan_backward_xy,  &
!                                  2*nxw, 2*nyw,      &
!                                  wwkk, wwxy,        &
!                                  FFTW_ESTIMATE )
!%%% setting for fftw3_omp %%%
      call dfftw_init_threads( ierr )
      call dfftw_plan_with_nthreads( omp_get_max_threads() )
      call dfftw_plan_dft_c2r_2d( plan_backward_xy,  &
                                  2*nxw, 2*nyw,      &
                                  wwkk, wwxy,        &
                                  FFTW_ESTIMATE )
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
      call dfftw_plan_with_nthreads( 1 )
      call dfftw_plan_dft_1d( plan_backward_x,   &
                              2*nxw,             &
                              wwkx, wwxx,        &
                              FFTW_BACKWARD,     &
                              FFTW_ESTIMATE )
  

  END SUBROUTINE fft_pre
  
  
!--------------------------------------
  SUBROUTINE fft_backward_xy ( gww, wwxy )
!--------------------------------------
!  Execution of FFT
  
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny)              :: gww
    real(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:2*nyw-1)             :: wwxy
  
    complex(kind=DP), dimension(0:nxw,0:2*nyw-1) :: wwkk
    integer :: mx, my
  
  
!$OMP parallel do
      do my = 0, 2*nyw-1
        do mx = 0, nxw
          wwkk(mx,my) = ( 0._DP, 0._DP )
        end do
      end do
!$OMP parallel do
      do my = 0, global_ny
        do mx = 0, nx
          wwkk(mx,my) = gww(mx,my)
        end do
      end do
!$OMP parallel do
      do my = 1, global_ny
        do mx = -nx, 0
          wwkk(-mx,2*nyw-my) = conjg( gww(mx,my) )
        end do
      end do
      mx = 0
        do my = 1, global_ny
          wwkk(mx,2*nyw-my) = conjg( wwkk(mx,my) )
        end do
  
      call dfftw_execute_dft_c2r( plan_backward_xy, wwkk, wwxy )

  
  END SUBROUTINE fft_backward_xy


!--------------------------------------
  SUBROUTINE fft_backward_x ( gww, wwxxky )
!--------------------------------------
!  Execution of FFT
  
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny)              :: gww
    complex(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:global_ny)           :: wwxxky
  
    complex(kind=DP), dimension(0:2*nxw-1) :: wwkx
    integer :: mx, my
  
  
!$OMP parallel do
      do my = 0, global_ny

        do mx = 0, nx
          wwkx(mx) = gww(mx,my)
        end do
        do mx = -nx, -1
          wwkx(2*nxw+mx) = gww(mx,my)
        end do
        do mx = nx+1, 2*nxw-nx-1
          wwkx(mx) = ( 0._DP, 0._DP )
        end do
  
        call dfftw_execute_dft( plan_backward_x, wwkx, wwxxky(:,my) )

      end do
  

  END SUBROUTINE fft_backward_x


END MODULE diag_fft


PROGRAM diag
!-------------------------------------------------------------------------------
!
!     Data diagnostics from binary output
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  use diag_header
  use diag_fft, only : fft_pre

  implicit none

  integer :: inum, loop
  character*3 :: cnum
  character*8 :: cloop

    call file_open

    if ( flag_cnt ) then
      loop = 0
      do inum = snum, enum
        if ( inum == snum ) then
          call init( inum )
          call fft_pre
        end if
        write( olog, * ) " #### Diagnostics of cnt. ####"
        write( cnum, fmt="(i3.3)" ) inum
        write( olog, * ) " ### inum = "//cnum//""
        write( cloop, fmt="(i8.8)" ) loop
        write( olog, * ) " ### loop_sta = "//cloop//""
        call diag_cnt( inum, loop )
        write( cloop, fmt="(i8.8)" ) loop-1
        write( olog, * ) " ### loop_end = "//cloop//""
        write( olog, * ) " #############################"
        write( olog, * ) ""
      end do
    end if

    if ( flag_mom ) then
      loop = 0
      do inum = snum, enum
        if ( inum == snum ) then
          call init( inum )
          call fft_pre
        end if
        write( olog, * ) " #### Diagnostics of mom. ####"
        write( cnum, fmt="(i3.3)" ) inum
        write( olog, * ) " ### inum = "//cnum//""
        write( cloop, fmt="(i8.8)" ) loop
        write( olog, * ) " ### loop_sta = "//cloop//""
        call diag_mom( inum, loop )
        write( cloop, fmt="(i8.8)" ) loop-1
        write( olog, * ) " ### loop_end = "//cloop//""
        write( olog, * ) " #############################"
        write( olog, * ) ""
      end do
    end if

    if ( flag_trn ) then
      loop = 0
      do inum = snum, enum
        if ( inum == snum ) then
          call init( inum )
          call fft_pre
        end if
        write( olog, * ) " #### Diagnostics of trn. ####"
        write( cnum, fmt="(i3.3)" ) inum
        write( olog, * ) " ### inum = "//cnum//""
        write( cloop, fmt="(i8.8)" ) loop
        write( olog, * ) " ### loop_sta = "//cloop//""
        call diag_trn( inum, loop )
        write( cloop, fmt="(i8.8)" ) loop-1
        write( olog, * ) " ### loop_end = "//cloop//""
        write( olog, * ) " #############################"
        write( olog, * ) ""
      end do
    end if

    call file_close

CONTAINS


!--------------------------------------
  SUBROUTINE file_open
!--------------------------------------
!     Open files
  
    use diag_header
  
    implicit none
  
      open( olog, file="./plt/log.dat" )
  
      if ( flag_cnt ) then

        open( pffinzv, file="./plot_ffinzv.gn" )
          write( pffinzv, * ) "#\!/usr/bin/gnuplot"
          write( pffinzv, * ) "# set terminal jpeg"
          write( pffinzv, * ) "  set contour base"
          write( pffinzv, * ) "  set cntrparam levels 20"
          write( pffinzv, * ) "  unset surface"
          write( pffinzv, * ) "  unset key"
          write( pffinzv, * ) "  set view 0, 0"
          write( pffinzv, * ) "  set xlabel 'Field-aligned coordinate z'"
          write( pffinzv, * ) "  set ylabel 'Parallel velocity vl'"

      end if

      if ( flag_mom ) then

        if ( flag_mominkxky ) then
          open( pmominkxky, file="./plot_mominkxky.gn" )
            write( pmominkxky, * ) "#\!/usr/bin/gnuplot"
            write( pmominkxky, * ) "# set terminal jpeg"
            write( pmominkxky, * ) "  set pm3d map"
            write( pmominkxky, * ) "  set size square"
            write( pmominkxky, * ) "  unset key"
            write( pmominkxky, * ) "  set xlabel 'Wave number kx'"
            write( pmominkxky, * ) "  set ylabel 'Wave number ky'"
            write( pmominkxky, * ) "  set palette rgbformulae 23,22,21  # cold"
            write( pmominkxky, * ) "# set palette rgbformulae 21,22,23  # hot"
            write( pmominkxky, * ) "# set palette rgbformulae 22,13,-31 # rainbow"
            write( pmominkxky, * ) "# set palette define ( -1 'blue', 0 'white', 1 'red' ) # positive and negative"
        end if
      
      end if

  
  END SUBROUTINE file_open
  
  
!--------------------------------------
  SUBROUTINE file_close
!--------------------------------------
!     Close files
  
    use diag_header
  
    implicit none
  
      close( olog )

      if ( flag_cnt ) then
        if ( flag_ffinzv    ) close( pffinzv )
      end if

      if ( flag_mom ) then
        if ( flag_mominkxky ) close( pmominkxky )
      end if
  
  
  END SUBROUTINE file_close
  
  
!--------------------------------------
  SUBROUTINE init( inum )
!--------------------------------------
!     Initial setting
  
    use diag_header
  
    implicit none
  
    integer, intent(in) :: inum
  
    character*3 :: cnum
    integer :: mx, my, iz, iv, im
  
    namelist /equib/ equib_type
    namelist /physp/ R0_Ln,  &    ! R0/Lns
                     R0_Lt,  &    ! R0/Lts
                        nu,  &    ! collision freq. 
                      Anum,  &    ! mass number
                      Znum,  &    ! charge number 
                       fcs,  &    ! charge-density fraction 
                       sgn,  &    ! signs of charge 
                       tau,  &    ! T-ratio Ts/T0, T0=reference temp.
                  lambda_i,  &    ! Debye/rho 
                      beta        ! mu0*ni*Ti/B^2

    namelist /nperi/ n_alp, n_tht
    namelist /confp/ r_minor, eps_r, eps_rnew,              &
                     q_0, q_d, s_hat,                       &
                     lprd, mprd, eps_hor, eps_mor, eps_por, &
                     rdeps00, rdeps1_0, rdeps1_10,          & 
                     rdeps2_10, rdeps3_10, malpha

  
      write( cnum, fmt="(i3.3)" ) inum
  
      open( 5, file="./gkvp_f0.30_namelist."//cnum, status="old", action="read" )
  
      read(5,nml=equib)
  
      read(5,nml=physp)
  
      read(5,nml=nperi)
  
      if( trim(equib_type) == "analytic"  .or.  &
          trim(equib_type) == "s-alpha"   .or.  &
          trim(equib_type) == "circ-MHD" ) then

        read(5,nml=confp)

        lprdm1   = lprd - 1.0_DP
        lprdp1   = lprd + 1.0_DP

        lmmq     = lprd   - mprd * q_0
        lmmqm1   = lprdm1 - mprd * q_0
        lmmqp1   = lprdp1 - mprd * q_0
  
      else
  
        write( olog, * ) " # Currently, this equilibrium is not available."
        stop
  
      end if
  
      lx       = r_minor * q_d / ( q_0 * s_hat )
      ly       = r_minor * pi  / ( q_0 * real( n_alp, kind=DP ) )
      lz       = real( n_tht, kind=DP ) * pi        ! total z-length
  
      kxmin    = pi / lx
      kymin    = pi / ly
  
      dz       = lz / real( global_nz, kind=DP )
  
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
        if( trim(equib_type) == "analytic" ) then
          gomg(iz) = 1._DP                                              &
                   - eps_r * ( cos( gzz(iz) )                           &
                           + eps_hor * cos( lmmq   * gzz(iz) - malpha ) &
                           + eps_mor * cos( lmmqm1 * gzz(iz) - malpha ) &
                           + eps_por * cos( lmmqp1 * gzz(iz) - malpha ) )
          grootg(iz) = q_0 / gomg(iz)
        else if( trim(equib_type) == "s-alpha"  ) then
        !  gomg(iz)   = 1._DP - eps_r * cos( gzz(iz) )
          gomg(iz)   = 1._DP / ( 1._DP + eps_r * cos( gzz(iz) ) ) ! s-alpha same as GENE
          grootg(iz) = q_0 / gomg(iz)
        else if( trim(equib_type) == "circ-MHD" ) then
          q_bar = dsqrt( 1._DP - eps_r**2 )*q_0
          theta = 2._DP*atan( sqrt( (1._DP+eps_r)/(1._DP-eps_r) ) &
                                             * tan(gzz(iz)/2._DP) )
          gomg(iz)   = sqrt( q_bar**2 + eps_r**2 ) &
                     / ( 1._DP + eps_r*cos( theta ) ) / q_bar
          grootg(iz) = (q_0**2/q_bar)*( 1._DP+eps_r*cos(theta) )**2
        end if
      end do
      cfsrf = 0._DP
      do iz = -global_nz, global_nz-1
        cfsrf = cfsrf + grootg(iz)
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
  
!$OMP parallel do
      do im = 0, global_nm
        do iv = 1, 2*global_nv
          do iz = -global_nz, global_nz-1
            gfmx(iz,iv,im) = exp( - 0.5_DP * gvl(iv)**2 - gomg(iz) * gmu(im) ) &
                           / sqrt( twopi**3 )
          end do
        end do
      end do
  
      close(5)
  
      write( olog, * ) "# DIAG parameters"
      write( olog, * ) "snum             =", snum            
      write( olog, * ) "enum             =", enum            
      write( olog, * ) "loopskip         =", loopskip        
      write( olog, * ) "flag_cnt         =", flag_cnt        
      write( olog, * ) "flag_ffinzv      =", flag_ffinzv     
      write( olog, * ) "flag_mom         =", flag_mom        
      write( olog, * ) "flag_mominkxky   =", flag_mominkxky  
      write( olog, * ) "flag_mominxy     =", flag_mominxy    
      write( olog, * ) "flag_mominz      =", flag_mominz     
      write( olog, * ) "flag_trninkxky   =", flag_trninkxky
      write( olog, * )
      write( olog, * ) "# GKV parameters"
      write( olog, * ) "nxw              =", nxw      
      write( olog, * ) "nyw              =", nyw      
      write( olog, * ) "nx               =", nx       
      write( olog, * ) "global_ny        =", global_ny
      write( olog, * ) "global_nz        =", global_nz
      write( olog, * ) "global_nv        =", global_nv
      write( olog, * ) "global_nm        =", global_nm
      write( olog, * ) "nzb              =", nzb      
      write( olog, * ) "nprocw           =", nprocw   
      write( olog, * ) "nprocz           =", nprocz   
      write( olog, * ) "nprocv           =", nprocv   
      write( olog, * ) "nprocm           =", nprocm   
      write( olog, * ) "nprocs           =", nprocs   
      write( olog, * ) "vmax             =", vmax     
      write( olog, * )

  
  END SUBROUTINE init
  

!--------------------------------------
  SUBROUTINE diag_cnt( inum, loop )
!--------------------------------------
!     Data diagnostics from *.cnt.*
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
    integer, intent(inout) :: loop
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir
    character*6 :: crank
    character*1 :: srank
    integer :: ios
  
    real(kind=DP) :: time
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: ff
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm,0:ns-1) :: gff
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor, &
               spc_rank, fft_rank, zsp_rank, vel_rank
    integer :: mx, my, gmy, iz, giz, iv, giv, im, gim
  

      write( cnum, fmt="(i3.3)" ) inum
  
  ! --- file open ---
      do ir = 0, nproc - 1

        rankg = ir
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
        scolor = mod( rankg, nprocz * nprocw )
        spc_rank = rankg / ( nprocz * nprocw )
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz ), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        fft_rank = mod( rank, nprocw )
        zsp_rank = mod( rank / nprocw, nprocz )
        vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
        
        write( crank, fmt="(i6.6)" ) ir
        write( srank, fmt="(i1.1)" ) ranks
        open( unit=100000+ir, file="./cnt/gkvp_f0.30."//crank//".cnt."//cnum,  &
              status="old", action="read", form="unformatted" )

      end do
  
  
  ! --- time-step loop ---
      do
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control
  
          do ir = 0, nproc - 1

            rankg = ir
            ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
            rank = mod( rankg, nproc / nprocs )
            scolor = mod( rankg, nprocz * nprocw )
            spc_rank = rankg / ( nprocz * nprocw )
            rankw = mod( rank,                      nprocw )
            rankz = mod( rank / nprocw,             nprocz )
            rankv = mod( rank / ( nprocw * nprocz ), nprocv )
            rankm = rank / ( nprocw * nprocz * nprocv )
            wcolor   = rank / nprocw
            zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
            vcolor   = rankz * nprocw  +  rankw
            fft_rank = mod( rank, nprocw )
            zsp_rank = mod( rank / nprocw, nprocz )
            vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
          
            read( unit=100000+ir, iostat=ios ) time, ff
            if ( ios /= 0 ) exit
!$OMP parallel do private(mx,my,iz,iv,im,gmy,giz,giv,gim)
            do im = 0, nm
              gim = (nm+1) * rankm + im
              do iv = 1, 2*nv
                giv = 2*nv * rankv + iv
                do iz = -nz, nz-1
                  giz = - global_nz + 2*nz * rankz + iz + nz
                  do my = 0, ny
                    gmy = ( ny+1 ) * rankw + my
                    if ( gmy <= global_ny ) then
                      do mx = -nx, nx
                        gff(mx,gmy,giz,giv,gim,ranks) = ff(mx,my,iz,iv,im)
                      end do
                    end if
                  end do
                end do
              end do
            end do
  
          end do
  
        else    !%%% loop control
  
          do ir = 0, nproc - 1
          
            read( unit=100000+ir, iostat=ios )
            if ( ios /= 0 ) exit
  
          end do
  
        end if  !%%% loop control
  
        if ( ios < 0 ) then
          write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                           ios, 100000+ir
          write( olog, * ) ""
          exit
        else if ( ios > 0 ) then
          write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                           ios, 100000+ir
          write( olog, * ) ""
          exit
        end if
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control

          mx = 0
          my = 6
          im = 4
          if ( flag_ffinzv      ) call ffinzv( mx, my, im, inum, loop, time, gff )

        end if  !%%% loop control
  
        loop = loop + 1
  
      end do
  
  
  ! --- file close ---
      do ir = 0, nproc - 1
        
        close( unit=100000+ir )
  
      end do
  
  
  END SUBROUTINE diag_cnt
  

!--------------------------------------
  SUBROUTINE ffinzv( mx, my, im, inum, loop, time, gff )
!--------------------------------------
!     Write |ff| in zz, vl space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: mx
    integer,          intent(in)                    :: my
    integer,          intent(in)                    :: im
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm,0:ns-1) :: gff
  
  ! --- local variables
  
    integer :: iconnect, connect_min, connect_max, mxw
    character*4 :: cmx, cmy, cim
    character*3 :: cnum
    character*8 :: cloop
    integer :: iz, iv
  
  
      write( cmx, fmt="(i4.4)" ) mx
      write( cmy, fmt="(i4.4)" ) my
      write( cim, fmt="(i4.4)" ) im
      write( cnum, fmt="(i3.3)" ) inum
      write( cloop, fmt="(i8.8)" ) loop
      open( offinzv, file="./plt/ffinzv_mx"//cmx//"my"//cmy//"im"//cim//"_t"//cloop//"."//cnum )
        write( offinzv, "(a17i17a17G17.7E3a17i17a17G17.7E3a17i17a17G17.7E3a17i17a17G17.7E3a17i17a17G17.7E3)" )  &
                                    "#           loop=",loop, "time=",time,           &
                                    "mx=",mx, "kx=",kx(mx), "my=",my, "ky=",gky(my),  &
                                    "im=",im, "mu=",gmu(im)
        write( offinzv, "(99a17)" ) "#               z","vl","|ff|"
  
        if ( dj(my) == 0 ) then
  
          do iv = 1, 2*global_nv
            do iz = -global_nz, global_nz-1
              write( offinzv, "(99G17.7E3)" ) gzz(iz), gvl(iv),  &
                                              abs( gff(mx,my,iz,iv,im,0:ns-1) )
            end do
            write( offinzv, * )
          end do
  
        else
  
          do iv = 1, 2*global_nv
            connect_min = int( ( nx + mx ) / abs( dj(my) ) )
            if ( connect_min .ne. 0 ) then
              do iconnect = connect_min, 1, -1
                mxw = mx+iconnect*dj(my)
                do iz = -global_nz, global_nz-1
                  write( offinzv, "(99G17.7E3)" )                    &
                       - twopi * real(iconnect) + gzz(iz), gvl(iv),  &
                       abs( ck(my)**iconnect * gff(mxw,my,iz,iv,im,0:ns-1) )
                end do
              end do
            end if
  
            connect_max = int( ( nx - mx ) / abs( dj(my) ) )
              do iconnect = 0, connect_max
                mxw = mx-iconnect*dj(my)
                do iz = -global_nz, global_nz-1
                  write( offinzv, "(99G17.7E3)" )                    &
                       + twopi * real(iconnect) + gzz(iz), gvl(iv),  &
                       abs( conjg( ck(my)**iconnect ) * gff(mxw,my,iz,iv,im,0:ns-1) )
                end do
              end do

            write( offinzv, * )
          end do
  
        end if

      close( offinzv )
  

  END SUBROUTINE ffinzv
  
  
!--------------------------------------
  SUBROUTINE diag_mom( inum, loop )
!--------------------------------------
!     Data diagnostics from *.phi.*, *.Al.*, *.mom.*
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
    integer, intent(inout) :: loop
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir
    character*6 :: crank
    character*1 :: srank
    integer :: ios
  
    real(kind=DP) :: time
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: dens, upara, ppara, pperp, qlpara, qlperp
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara, &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor, &
               spc_rank, fft_rank, zsp_rank, vel_rank
    integer :: mx, my, gmy, iz, giz
  
      write( cnum, fmt="(i3.3)" ) inum
  
  ! --- file open ---
      do ir = 0, nproc - 1

        rankg = ir
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
        scolor = mod( rankg, nprocz * nprocw )
        spc_rank = rankg / ( nprocz * nprocw )
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz ), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        fft_rank = mod( rank, nprocw )
        zsp_rank = mod( rank / nprocw, nprocz )
        vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
        
        write( crank, fmt="(i6.6)" ) ir
        write( srank, fmt="(i1.1)" ) ranks

        if ( ranks == 0 .and. vel_rank == 0 ) then
          open( unit=200000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".phi."//cnum,  &
                status="old", action="read", form="unformatted" )
          open( unit=300000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".Al."//cnum,  &
                status="old", action="read", form="unformatted" )
        end if

        if ( vel_rank == 0 ) then
          open( unit=400000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".mom."//cnum,  &
                status="old", action="read", form="unformatted" )
        end if

      end do
  
  
  ! --- time-step loop ---
      do
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control
  
          do ir = 0, nproc - 1
  
            rankg = ir
            ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
            rank = mod( rankg, nproc / nprocs )
            scolor = mod( rankg, nprocz * nprocw )
            spc_rank = rankg / ( nprocz * nprocw )
            rankw = mod( rank,                      nprocw )
            rankz = mod( rank / nprocw,             nprocz )
            rankv = mod( rank / ( nprocw * nprocz ), nprocv )
            rankm = rank / ( nprocw * nprocz * nprocv )
            wcolor   = rank / nprocw
            zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
            vcolor   = rankz * nprocw  +  rankw
            fft_rank = mod( rank, nprocw )
            zsp_rank = mod( rank / nprocw, nprocz )
            vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
            
            if ( ranks == 0 .and. vel_rank == 0 ) then
              read( unit=200000+ir, iostat=ios ) time, phi
              read( unit=300000+ir, iostat=ios ) time, Al
              if ( ios /= 0 ) exit
!$OMP parallel do private(mx,my,iz,gmy,giz)
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
  
            if ( vel_rank == 0 ) then
              read( unit=400000+ir, iostat=ios ) time, dens, upara, ppara, pperp, qlpara, qlperp
              if ( ios /= 0 ) exit
!$OMP parallel do private(mx,my,iz,gmy,giz)
              do iz = -nz, nz-1
                giz = - global_nz + 2*nz * rankz + iz + nz
                do my = 0, ny
                  gmy = ( ny+1 ) * rankw + my
                  if ( gmy <= global_ny ) then
                    do mx = -nx, nx
                        gdens(mx,gmy,giz,ranks) =   dens(mx,my,iz)
                       gupara(mx,gmy,giz,ranks) =  upara(mx,my,iz)
                       gppara(mx,gmy,giz,ranks) =  ppara(mx,my,iz)
                       gpperp(mx,gmy,giz,ranks) =  pperp(mx,my,iz)
                      gqlpara(mx,gmy,giz,ranks) = qlpara(mx,my,iz)
                      gqlperp(mx,gmy,giz,ranks) = qlperp(mx,my,iz)
                    end do
                  end if
                end do
              end do
            end if
  
          end do
  
        else    !%%% loop control
  
          do ir = 0, nproc - 1

            rankg = ir
            ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
            rank = mod( rankg, nproc / nprocs )
            scolor = mod( rankg, nprocz * nprocw )
            spc_rank = rankg / ( nprocz * nprocw )
            rankw = mod( rank,                      nprocw )
            rankz = mod( rank / nprocw,             nprocz )
            rankv = mod( rank / ( nprocw * nprocz ), nprocv )
            rankm = rank / ( nprocw * nprocz * nprocv )
            wcolor   = rank / nprocw
            zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
            vcolor   = rankz * nprocw  +  rankw
            fft_rank = mod( rank, nprocw )
            zsp_rank = mod( rank / nprocw, nprocz )
            vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
            
            if ( ranks == 0 .and. vel_rank == 0 ) then
              read( unit=200000+ir, iostat=ios )
              if ( ios /= 0 ) exit
              read( unit=300000+ir, iostat=ios )
              if ( ios /= 0 ) exit
            end if
  
            if ( vel_rank == 0 ) then
              read( unit=400000+ir, iostat=ios )
              if ( ios /= 0 ) exit
            end if
  
          end do
  
        end if  !%%% loop control
  
        if ( ios < 0 ) then
          write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                           ios, 200000+ir, 300000+ir, 400000+ir
          write( olog, * ) ""
          exit
        else if ( ios > 0 ) then
          write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                           ios, 200000+ir, 300000+ir, 400000+ir
          write( olog, * ) ""
          exit
        end if
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control

          if ( flag_mominkxky ) then
            call mominkxky( inum, loop, time,  &
                            gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
          end if
          if ( flag_mominxy ) then
            iz = 0
            call mominxy( iz, inum, loop, time,  &
                          gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
          end if
          if ( flag_mominz ) then
            mx = 0
            my = 6
            call mominz( mx, my, inum, loop, time,  &
                         gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
          end if

          if ( mod(loop,10*loopskip)==0 ) write(*,*) "inum, loop, time =",inum,loop,time

        end if  !%%% loop control
  
        loop = loop + 1
  
      end do
  
  
  ! --- file close ---
      do ir = 0, nproc - 1

        rankg = ir
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
        scolor = mod( rankg, nprocz * nprocw )
        spc_rank = rankg / ( nprocz * nprocw )
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz ), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        fft_rank = mod( rank, nprocw )
        zsp_rank = mod( rank / nprocw, nprocz )
        vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
        
        if ( ranks == 0 .and. vel_rank == 0 ) then
          close( unit=200000+ir )
          close( unit=300000+ir )
        end if
  
        if ( vel_rank == 0 ) then
          close( unit=400000+ir )
        end if
  
      end do
  
  
  END SUBROUTINE diag_mom
  
  
!--------------------------------------
  SUBROUTINE mominkxky( inum, loop, time, gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!--------------------------------------
!     Write phi in kx, ky space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara, &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
  
  ! --- local variables
  
    real(kind=DP), dimension(-nx:nx,0:global_ny)           :: phi2, Al2
    real(kind=DP), dimension(-nx:nx,0:global_ny,0:ns-1)    :: dens2, upara2, pres2, qpara2,  &
                                                              phidens, uparaAl, phipres, qparaAl
    real(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wr3
    character*3 :: cnum
    character*8 :: cloop
    integer :: mx, my, iz, is
  
  
!$OMP parallel do
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            wr3(mx,my,iz) = real( gphi(mx,my,iz) * conjg( gphi(mx,my,iz) ), kind=DP )
          end do
        end do
      end do
      call thet_ave_r ( wr3, phi2 )
!$OMP parallel do
      do iz = -global_nz, global_nz-1
        do my = 0, global_ny
          do mx = -nx, nx
            wr3(mx,my,iz) = real( gAl(mx,my,iz) * conjg( gAl(mx,my,iz) ), kind=DP )
          end do
        end do
      end do
      call thet_ave_r ( wr3, Al2 )
      do is = 0, ns-1
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( gdens(mx,my,iz,is) * conjg( gdens(mx,my,iz,is) ), kind=DP )
            end do
          end do
        end do
        call thet_ave_r ( wr3, dens2(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( gupara(mx,my,iz,is) * conjg( gupara(mx,my,iz,is) ), kind=DP )
            end do
          end do
        end do
        call thet_ave_r ( wr3, upara2(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( ( gppara(mx,my,iz,is) + gpperp(mx,my,iz,is) )  &
                                    * conjg( gppara(mx,my,iz,is) + gpperp(mx,my,iz,is) ), kind=DP )
            end do
          end do
        end do
        call thet_ave_r ( wr3, pres2(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = real( ( gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is) )  &
                                    * conjg( gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is) ), kind=DP )
            end do
          end do
        end do
        call thet_ave_r ( wr3, qpara2(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = aimag( gphi(mx,my,iz) * conjg( gdens(mx,my,iz,is) ) )
            end do
          end do
        end do
        call thet_ave_r ( wr3, phidens(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = aimag( gupara(mx,my,iz,is) * conjg( gAl(mx,my,iz) ) )
            end do
          end do
        end do
        call thet_ave_r ( wr3, uparaAl(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = aimag( gphi(mx,my,iz)  &
                                     * conjg( gppara(mx,my,iz,is) + gpperp(mx,my,iz,is) ) )
            end do
          end do
        end do
        call thet_ave_r ( wr3, phipres(:,:,is) )
!$OMP parallel do
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wr3(mx,my,iz) = aimag( ( gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is) )  &
                                     * conjg( gAl(mx,my,iz) ) )
            end do
          end do
        end do
        call thet_ave_r ( wr3, qparaAl(:,:,is) )
      end do

      write( cnum, fmt="(i3.3)" ) inum
      write( cloop, fmt="(i8.8)" ) loop
      open( omominkxky, file="./plt/mominkxky_t"//cloop//"."//cnum )
        write( omominkxky, "(a17i17a17G17.7E3)" )  &
                                       "#           loop=",loop, "time=",time
        write( omominkxky, "(99a17)" ) "#              kx","ky","||phi||","||Al||",    &
                                       "||dens||","||upara||","||pres||","||qpara||",  &
                                       "Im[<phi,dens>]","Im[<upara,Al>]",              &
                                       "Im[<phi,pres>]","Im[<qpara,Al>]"
        do my = 0, global_ny
          do mx = -nx, nx
            write( omominkxky, "(99G17.7E3)" ) kx(mx), gky(my), phi2(mx,my), Al2(mx,my),  &
                                             (   dens2(mx,my,is),  &
                                                upara2(mx,my,is),  &
                                                 pres2(mx,my,is),  &
                                                qpara2(mx,my,is),  &
                                               phidens(mx,my,is),  &
                                               uparaAl(mx,my,is),  &
                                               phipres(mx,my,is),  &
                                               qparaAl(mx,my,is), is=0, ns-1 )
          end do
          write( omominkxky, * )
        end do
      close( omominkxky )
  
      write( pmominkxky, * ) "# set output './plt/mominkxky_t"//cloop//"."//cnum//".jpg'"
      write( pmominkxky, * ) "  set title 'time =",time,"'"
      write( pmominkxky, * ) "  splot './plt/mominkxky_t"//cloop//"."//cnum//"' u 1:2:3"
      write( pmominkxky, * ) "  pause 0.2"

  
  END SUBROUTINE mominkxky
  
  
!--------------------------------------
  SUBROUTINE mominxy( iz, inum, loop, time, gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!--------------------------------------
!     Write phi in x, y space
  
    use diag_header
    use diag_fft, only : fft_backward_xy
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: iz
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara, &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
  
  ! --- local variables
  
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1)        :: phixy, Alxy
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1,0:ns-1) :: densxy, uparaxy, presxy, qparaxy
    complex(kind=DP), dimension(-nx:nx,0:global_ny) :: wc2
    character*4 :: ciz
    character*3 :: cnum
    character*8 :: cloop
    integer :: mx, my, is
  
  
!$OMP parallel do
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          phixy(mx,my) = 0._DP
           Alxy(mx,my) = 0._DP
        end do
      end do
      do is = 0, ns-1
!$OMP parallel do
        do my = 0, 2*nyw-1
          do mx = 0, 2*nxw-1
             densxy(mx,my,is) = 0._DP
            uparaxy(mx,my,is) = 0._DP
             presxy(mx,my,is) = 0._DP
            qparaxy(mx,my,is) = 0._DP
          end do
        end do
      end do

      call fft_backward_xy( gphi(:,:,iz), phixy )
      call fft_backward_xy(  gAl(:,:,iz),  Alxy )
      do is = 0, ns-1
        call fft_backward_xy(  gdens(:,:,iz,is),  densxy(:,:,is) )
        call fft_backward_xy( gupara(:,:,iz,is), uparaxy(:,:,is) )
!$OMP parallel do
        do my = 0, global_ny
          do mx = -nx, nx
            wc2(mx,my) = gppara(mx,my,iz,is) + gpperp(mx,my,iz,is)
          end do
        end do
        call fft_backward_xy( wc2, presxy(:,:,is) )
!$OMP parallel do
        do my = 0, global_ny
          do mx = -nx, nx
            wc2(mx,my) = gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is)
          end do
        end do
        call fft_backward_xy( wc2, qparaxy(:,:,is) )
      end do

      write( ciz, fmt="(i4.4)" ) iz
      write( cnum, fmt="(i3.3)" ) inum
      write( cloop, fmt="(i8.8)" ) loop
      open( omominxy, file="./plt/mominxy_z"//ciz//"_t"//cloop//"."//cnum )
        write( omominxy, "(a17i17a17G17.7E3a17i17a17G17.7E3)" )  &
                                     "#           loop=",loop, "time=",time, "iz=",iz, "zz=",gzz(iz)
        write( omominxy, "(99a17)" ) "#               x","y","phi","Al",  &
                                     "dens","upara","pres","qpara"
        do my = 0, 2*nyw-1
          do mx = 0, 2*nxw-1
            write( omominxy, "(99G17.7E3)" ) xx(mx), yy(my), phixy(mx,my), Alxy(mx,my),  &
                                           (  densxy(mx,my,is),  &
                                             uparaxy(mx,my,is),  &
                                              presxy(mx,my,is),  &
                                             qparaxy(mx,my,is), is=0, ns-1 )
          end do
          write( omominxy, * )
        end do
      close( omominxy )
 
  
  END SUBROUTINE mominxy
  
  
!--------------------------------------
  SUBROUTINE mominz( mx, my, inum, loop, time, gphi, gAl, gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!--------------------------------------
!     Write phi in zz space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: mx
    integer,          intent(in)                    :: my
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara, &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
  
  
  ! --- local variables
  
    integer :: iconnect, connect_min, connect_max, mxw
    complex(kind=DP) :: wc
    character*4 :: cmx, cmy
    character*3 :: cnum
    character*8 :: cloop
    integer :: iz, is
  
  
      write( cmx, fmt="(i4.4)" ) mx
      write( cmy, fmt="(i4.4)" ) my
      write( cnum, fmt="(i3.3)" ) inum
      write( cloop, fmt="(i8.8)" ) loop
      open( omominz, file="./plt/mominz_mx"//cmx//"my"//cmy//"_t"//cloop//"."//cnum )
        write( omominz, "(a17i17a17G17.7E3a17i17a17G17.7E3a17i17a17G17.7E3)" )  &
                                    "#           loop=",loop, "time=",time,  &
                                    "mx=",mx, "kx=",kx(mx), "my=",my, "ky=",gky(my)
        write( omominz, "(99a17)" ) "#               z",  &
                                    "Re[phi]","Im[phi]","Re[Al]","Im[Al]",  &
                                    "Re[dens]","Im[dens]","Re[upara]","Im[upara]",  &
                                    "Re[pres]","Im[pres]","Re[qpara]","Im[qpara]"

        if ( dj(my) == 0 ) then
  
            do iz = -global_nz, global_nz-1
              write( omominz, "(99G17.7E3)" ) gzz(iz),  &
                                              gphi(mx,my,iz),  &
                                              gAl(mx,my,iz),   &
                                            ( gdens(mx,my,iz,is),                         &
                                              gupara(mx,my,iz,is),                        &
                                              gppara(mx,my,iz,is)+gpperp(mx,my,iz,is),    &
                                              gqlpara(mx,my,iz,is)+gqlperp(mx,my,iz,is),  &
                                                is=0, ns-1 )
            end do
  
        else
  
          connect_min = int( ( nx + mx ) / abs( dj(my) ) )
          if ( connect_min .ne. 0 ) then
            do iconnect = connect_min, 1, -1
              mxw = mx+iconnect*dj(my)
              wc = ck(my)**iconnect
              do iz = -global_nz, global_nz-1
                write( omominz, "(99G17.7E3)" ) - twopi * real(iconnect) + gzz(iz),  &
                                                wc*gphi(mxw,my,iz),  &
                                                wc*gAl(mxw,my,iz),   &
                                              ( wc*gdens(mxw,my,iz,is),                            &
                                                wc*gupara(mxw,my,iz,is),                           &
                                                wc*(gppara(mxw,my,iz,is)+gpperp(mxw,my,iz,is)),    &
                                                wc*(gqlpara(mxw,my,iz,is)+gqlperp(mxw,my,iz,is)),  &
                                                  is=0, ns-1 )
              end do
            end do
          end if
  
          connect_max = int( ( nx - mx ) / abs( dj(my) ) )
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(my)
              wc = conjg( ck(my)**iconnect )
              do iz = -global_nz, global_nz-1
                write( omominz, "(99G17.7E3)" ) + twopi * real(iconnect) + gzz(iz),  &
                                                wc*gphi(mxw,my,iz),  &
                                                wc*gAl(mxw,my,iz),   &
                                              ( wc*gdens(mxw,my,iz,is),                            &
                                                wc*gupara(mxw,my,iz,is),                           &
                                                wc*(gppara(mxw,my,iz,is)+gpperp(mxw,my,iz,is)),    &
                                                wc*(gqlpara(mxw,my,iz,is)+gqlperp(mxw,my,iz,is)),  &
                                                  is=0, ns-1 )
              end do
            end do
 
        end if

        write( omominz, * )

      close( omominz )
  
  
  END SUBROUTINE mominz
  

!--------------------------------------
  SUBROUTINE diag_trn( inum, loop )
!--------------------------------------
!     Data diagnostics from *.trn.*
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
    integer, intent(inout) :: loop
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir
    character*6 :: crank
    character*1 :: srank
    integer :: ios
  
    real(kind=DP) :: time
    real(kind=DP),  &
      dimension(-nx:nx,0:ny)          :: entrpy, neint, nmint, dcd
    real(kind=DP),  &
      dimension(-nx:nx,0:global_ny,0:ns-1)      :: gentrpy, gneint, gnmint, gdcd
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor, &
               spc_rank, fft_rank, zsp_rank, vel_rank
    integer :: mx, my, gmy, iz, giz
  
      write( cnum, fmt="(i3.3)" ) inum
  
  ! --- file open ---
      do ir = 0, nproc - 1

        rankg = ir
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
        scolor = mod( rankg, nprocz * nprocw )
        spc_rank = rankg / ( nprocz * nprocw )
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz ), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        fft_rank = mod( rank, nprocw )
        zsp_rank = mod( rank / nprocw, nprocz )
        vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
        
        write( crank, fmt="(i6.6)" ) ir
        write( srank, fmt="(i1.1)" ) ranks

        if ( zsp_rank == 0 .and. vel_rank == 0 ) then
          open( unit=500000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".trn."//cnum,  &
                status="old", action="read", form="unformatted" )
        end if
  
      end do
  
  
  ! --- time-step loop ---
      do
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control
  
          do ir = 0, nproc - 1
  
            rankg = ir
            ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
            rank = mod( rankg, nproc / nprocs )
            scolor = mod( rankg, nprocz * nprocw )
            spc_rank = rankg / ( nprocz * nprocw )
            rankw = mod( rank,                      nprocw )
            rankz = mod( rank / nprocw,             nprocz )
            rankv = mod( rank / ( nprocw * nprocz ), nprocv )
            rankm = rank / ( nprocw * nprocz * nprocv )
            wcolor   = rank / nprocw
            zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
            vcolor   = rankz * nprocw  +  rankw
            fft_rank = mod( rank, nprocw )
            zsp_rank = mod( rank / nprocw, nprocz )
            vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
            
            if ( zsp_rank == 0 .and. vel_rank == 0 ) then
              read( unit=500000+ir, iostat=ios ) time, entrpy, neint, nmint, dcd
              if ( ios /= 0 ) exit
!$OMP parallel do private(mx,my,gmy)
              do my = 0, ny
                gmy = ( ny+1 ) * rankw + my
                if ( gmy <= global_ny ) then
                  do mx = -nx, nx
                    gentrpy(mx,gmy,ranks) = entrpy(mx,my)
                     gneint(mx,gmy,ranks) =  neint(mx,my)
                     gnmint(mx,gmy,ranks) =  nmint(mx,my)
                       gdcd(mx,gmy,ranks) =    dcd(mx,my)
                  end do
                end if
              end do
            end if
  
          end do
  
        else    !%%% loop control
  
          do ir = 0, nproc - 1

            rankg = ir
            ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
            rank = mod( rankg, nproc / nprocs )
            scolor = mod( rankg, nprocz * nprocw )
            spc_rank = rankg / ( nprocz * nprocw )
            rankw = mod( rank,                      nprocw )
            rankz = mod( rank / nprocw,             nprocz )
            rankv = mod( rank / ( nprocw * nprocz ), nprocv )
            rankm = rank / ( nprocw * nprocz * nprocv )
            wcolor   = rank / nprocw
            zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
            vcolor   = rankz * nprocw  +  rankw
            fft_rank = mod( rank, nprocw )
            zsp_rank = mod( rank / nprocw, nprocz )
            vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
            
            if ( zsp_rank == 0 .and. vel_rank == 0 ) then
              read( unit=500000+ir, iostat=ios )
              if ( ios /= 0 ) exit
            end if
    
          end do
  
        end if  !%%% loop control
  
        if ( ios < 0 ) then
          write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                           ios, 500000+ir
          write( olog, * ) ""
          exit
        else if ( ios > 0 ) then
          write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                           ios, 500000+ir
          write( olog, * ) ""
          exit
        end if
  
        if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then  !%%% loop control

          if ( flag_trninkxky ) then
            call trninkxky( inum, loop, time,  &
                              gentrpy, gneint, gnmint, gdcd )
          end if

          if ( mod(loop,10*loopskip)==0 ) write(*,*) "inum, loop, time =",inum,loop,time

        end if  !%%% loop control
  
        loop = loop + 1
  
      end do
  
  
  ! --- file close ---
      do ir = 0, nproc - 1

        rankg = ir
        ranks = rankg / ( nprocz * nprocv * nprocm * nprocw )
        rank = mod( rankg, nproc / nprocs )
        scolor = mod( rankg, nprocz * nprocw )
        spc_rank = rankg / ( nprocz * nprocw )
        rankw = mod( rank,                      nprocw )
        rankz = mod( rank / nprocw,             nprocz )
        rankv = mod( rank / ( nprocw * nprocz ), nprocv )
        rankm = rank / ( nprocw * nprocz * nprocv )
        wcolor   = rank / nprocw
        zcolor   = ( rank / ( nprocw * nprocz ) ) * nprocw + rankw
        vcolor   = rankz * nprocw  +  rankw
        fft_rank = mod( rank, nprocw )
        zsp_rank = mod( rank / nprocw, nprocz )
        vel_rank = mod( rank / ( nprocw * nprocz ), nprocv * nprocm )
        
        if ( zsp_rank == 0 .and. vel_rank == 0 ) then
          close( unit=500000+ir )
        end if
  
      end do
  
  
  END SUBROUTINE diag_trn

  
!--------------------------------------
  SUBROUTINE trninkxky( inum, loop, time, gentrpy, gneint, gnmint, gdcd )
!--------------------------------------
!     Write nonlinear transfer in kx, ky

    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)  :: inum
    integer,          intent(in)  :: loop
    real(kind=DP),    intent(in)  :: time
    real(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,0:ns-1)      :: gentrpy, gneint, gnmint, gdcd
  
  ! --- local variables
  
    character*3 :: cnum
    character*8 :: cloop
    integer :: mx, my, is
  
  
      write( cnum, fmt="(i3.3)" ) inum
      write( cloop, fmt="(i8.8)" ) loop
      open( otrninkxky, file="./plt/trninkxky_t"//cloop//"."//cnum )
        write( otrninkxky, "(a17i17a17G17.7E3)" )  &
                                         "#           loop=",loop, "time=",time
        write( otrninkxky, "(99a17)" ) "#              kx","ky","entrpy","neint","nmint","dcd"
        do my = 0, global_ny
          do mx = -nx, nx
            write( otrninkxky, "(99G17.7E3)" ) kx(mx), gky(my),    &
                                             ( gentrpy(mx,my,is),  &
                                                gneint(mx,my,is),  &
                                                gnmint(mx,my,is),  &
                                                  gdcd(mx,my,is), is=0, ns-1 )
          end do
          write( otrninkxky, * )
        end do

      close( otrninkxky )
  
  
  END SUBROUTINE trninkxky


!--------------------------------------
  SUBROUTINE thet_ave_z ( wn, wa )
!--------------------------------------
!     average of a complex variable wn in the theta space

    use diag_header
  
    implicit none
  
  ! --- arguments
  
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wn
  
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny)                        :: wa
  
  ! --- local variables
  
    real(kind=DP) :: fct
    integer ::  mx, my, iz
  
  
      wa   = ( 0._DP, 0._DP )
  
!$OMP parallel do private(fct) reduction(+:wa)
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / cfsrf
        do my = 0, global_ny
          do mx = -nx, nx
            wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
          end do
        end do
      end do
  
  
  END SUBROUTINE thet_ave_z


!--------------------------------------
  SUBROUTINE thet_ave_r ( wn, wa )
!--------------------------------------
!     average of a real variable wn in the theta space

    use diag_header
  
    implicit none
  
  ! --- arguments
  
    real(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wn
  
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny)                        :: wa
  
  ! --- local variables
  
    real(kind=DP) :: fct
    integer ::  mx, my, iz
  
  
      wa   = 0._DP
  
!$OMP parallel do private(fct) reduction(+:wa)
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / cfsrf
        do my = 0, global_ny
          do mx = -nx, nx
            wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
          end do
        end do
      end do
  
  
  END SUBROUTINE thet_ave_r


END PROGRAM diag
