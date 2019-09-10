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
  logical, parameter :: flag_ffinzv      = .false.
  logical, parameter :: flag_mom = .true.
  logical, parameter :: flag_phiinkxky   = .false.
  logical, parameter :: flag_phiinxy     = .true.
  logical, parameter :: flag_phiinz      = .false.
  logical, parameter :: flag_Alinkxky    = .false.
  logical, parameter :: flag_Alinxy      = .false.
  logical, parameter :: flag_Alinz       = .false.
  logical, parameter :: flag_fluxinkxky  = .false.
  logical, parameter :: flag_fluxinky    = .false.
  logical, parameter :: flag_transinkxky = .false.
!%%%%%%%%%%%%%%%%%%%%%%

!%%% GKV parameters %%%
  integer, parameter :: nxw = 512, nyw = 512
  integer, parameter :: nx = 320, global_ny = 319  ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 32, global_nv = 48, global_nm = 15
  integer, parameter :: nzb = 2
  integer, parameter :: nprocw = 32, nprocz = 8, nprocv = 12, nprocm = 2, nprocs = 2
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
                     offinzv = 20, &
                  ophiinkxky = 40, &
                    ophiinxy = 41, &
                     ophiinz = 42, &
                   oAlinkxky = 50, &
                     oAlinxy = 51, &
                      oAlinz = 52, &
                 ofluxinkxky = 60, &
                   ofluxinky = 61, &
                otransinkxky = 70


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


      call dfftw_plan_dft_c2r_2d( plan_backward_xy,  &
                                  2*nxw, 2*nyw,      &
                                  wwkk, wwxy,        &
                                  FFTW_ESTIMATE )
  
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
  
  
      do my = 0, 2*nyw-1
        do mx = 0, nxw
          wwkk(mx,my) = ( 0._DP, 0._DP )
        end do
      end do
      do my = 0, global_ny
        do mx = 0, nx
          wwkk(mx,my) = gww(mx,my)
        end do
      end do
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

      if ( flag_cnt ) call diag_cnt( inum )
      if ( flag_mom ) call diag_mom( inum )

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
  
    use diag_header
  
    implicit none
  
      open( olog, file="./plt/log.dat" )
  
      if ( flag_cnt ) then

        if ( flag_ffinzv ) then
          open( offinzv, file="./plt/ffinzvtime.dat" )
          write( offinzv, "(99a17)" ) "#              zz","vl","time","|ff|"
        end if

      end if

      if ( flag_mom ) then

        if ( flag_phiinkxky ) then
          open( ophiinkxky, file="./plt/phiinkxkytime.dat" )
          write( ophiinkxky, "(99a17)" ) "#              kx","ky","time","|phi|"
        end if
      
        if ( flag_phiinxy ) then
          open( ophiinxy, file="./plt/phiinxytime.dat" )
          write( ophiinxy, "(99a17)" ) "#               x","y","time","phi"
        end if
      
        if ( flag_phiinz ) then
          open( ophiinz, file="./plt/phiinztime.dat" )
          write( ophiinz, "(99a17)" ) "#               z","time","Re[phi]","Im[phi]"
        end if
      
        if ( flag_Alinkxky ) then
          open( oAlinkxky, file="./plt/Alinkxkytime.dat" )
          write( oAlinkxky, "(99a17)" ) "#              kx","ky","time","|Al|"
        end if
      
        if ( flag_Alinxy ) then
          open( oAlinxy, file="./plt/Alinxytime.dat" )
          write( oAlinxy, "(99a17)" ) "#               x","y","time","Al"
        end if
      
        if ( flag_Alinz ) then
          open( oAlinz, file="./plt/Alinztime.dat" )
          write( oAlinz, "(99a17)" ) "#               z","time","Re[Al]","Im[Al]"
        end if
      
        if ( flag_fluxinkxky ) then
          open( ofluxinkxky, file="./plt/fluxinkxkytime.dat" )
          write( ofluxinkxky, "(99a17)" ) "#              kx","ky","time", &
                                          "p_flux_es(ns)","p_flux_em(ns)", &
                                          "h_flux_es(ns)","h_flux_em(ns)"
        end if
    
        if ( flag_fluxinky ) then
          open( ofluxinky, file="./plt/fluxinkytime.dat" )
            write( ofluxinky, "(99a17)" ) "#              ky", "time",  &
                                     "p_flux_es(ns)", "p_flux_em(ns)",  &
                                     "h_flux_es(ns)", "h_flux_em(ns)",  &
                                     "p_ampl_es(ns)", "p_ampl_em(ns)",  &
                                     "h_ampl_es(ns)", "h_ampl_em(ns)"
        end if
    
        if ( flag_transinkxky ) then
          open( otransinkxky, file="./plt/transinkxkytime.dat" )
          write( otransinkxky, "(99a17)" ) "#              kx","ky","time", &
                                           "entrpy(ns)","neint(ns)","nmint(ns)","dcd(ns)"
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
        if ( flag_ffinzv      ) close( offinzv )
      end if

      if ( flag_mom ) then
        if ( flag_phiinkxky   ) close( ophiinkxky )
        if ( flag_phiinxy     ) close( ophiinxy )
        if ( flag_phiinz      ) close( ophiinz )
        if ( flag_Alinkxky    ) close( oAlinkxky )
        if ( flag_Alinxy      ) close( oAlinxy )
        if ( flag_Alinz       ) close( oAlinz )
        if ( flag_fluxinkxky  ) close( ofluxinkxky )
        if ( flag_fluxinky    ) close( ofluxinky )
        if ( flag_transinkxky ) close( otransinkxky )
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
  
      if( trim(equib_type) == "analytic"  .OR.  &
          trim(equib_type) == "s-alpha"   .OR.  &
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
  SUBROUTINE diag_cnt( inum )
!--------------------------------------
!     Data diagnostics from *.cnt.*
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir
    character*6 :: crank
    character*1 :: srank
    integer :: loop, ios
  
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
      write( olog, * ) " # Diagnostics of cnt."
  
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
      !  write( olog, * ) " # Opened file and the unit-ID are ",  &
      !                   "./cnt/gkvp_f0.30."//crank//".cnt."//cnum, 100000+ir

      end do
  
  
  ! --- time-step loop ---
      loop = 0
      do
  
    if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then            !%%% to restart
  
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
  
    else                             !%%% to restart
  
        do ir = 0, nproc - 1
          
          read( unit=100000+ir, iostat=ios )
          if ( ios /= 0 ) exit
  
        end do
  
    end if                           !%%% to restart
  
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
  
    if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then            !%%% to restart

        if ( flag_ffinzv      ) call ffinzv( inum, loop, time, gff )
      !  write( olog, * ) " # Data output at loop, time = ", loop, time

    end if                           !%%% to restart
  
        loop = loop + 1
  
      end do
  
  
  ! --- file close ---
      do ir = 0, nproc - 1
        
        close( unit=100000+ir )
  
      end do
  
  
  END SUBROUTINE diag_cnt
  

!--------------------------------------
  SUBROUTINE ffinzv( inum, loop, time, gff )
!--------------------------------------
!     Write |ff| in zz, vl space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,1:2*global_nv,0:global_nm,0:ns-1) :: gff
  
  ! --- local variables
  
    integer :: iconnect, connect_min, connect_max, mxw
    integer :: mx, my, iz, iv, im
  
  
      mx = 0
      my = 6
      im = 4
  
      if ( dj(my) == 0 ) then
  
        do iv = 1, 2*global_nv
          do iz = -global_nz, global_nz-1
            write( offinzv, "(99G17.7E3)" ) gzz(iz), gvl(iv), time,  &
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
                write( offinzv, "(99G17.7E3)" )                          &
                     - twopi * real(iconnect) + gzz(iz), gvl(iv), time,  &
                     abs( ck(my)**iconnect * gff(mxw,my,iz,iv,im,0:ns-1) )
              end do
            end do
          end if
  
          connect_max = int( ( nx - mx ) / abs( dj(my) ) )
            do iconnect = 0, connect_max
              mxw = mx-iconnect*dj(my)
              do iz = -global_nz, global_nz-1
                write( offinzv, "(99G17.7E3)" )                          &
                     + twopi * real(iconnect) + gzz(iz), gvl(iv), time,  &
                     abs( conjg( ck(my)**iconnect ) * gff(mxw,my,iz,iv,im,0:ns-1) )
              end do
            end do

          write( offinzv, * )
        end do
  
      end if

  
  END SUBROUTINE ffinzv
  
  
!--------------------------------------
  SUBROUTINE diag_mom( inum )
!--------------------------------------
!     Data diagnostics from *.phi.*, *.Al.*, *.mom.*
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer, intent(in) :: inum
  
  ! --- local variables
  
    character*3 :: cnum
    integer :: ir
    character*6 :: crank
    character*1 :: srank
    integer :: loop, ios
  
    real(kind=DP) :: time
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: phi, Al
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: dens, upara, ppara, pperp, qlpara, qlperp
    real(kind=DP),  &
      dimension(-nx:nx,0:ny)          :: entrpy, neint, nmint, dcd
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi, gAl
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara, &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
    real(kind=DP),  &
      dimension(-nx:nx,0:global_ny,0:ns-1)      :: gentrpy, gneint, gnmint, gdcd
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor, &
               spc_rank, fft_rank, zsp_rank, vel_rank
    integer :: mx, my, gmy, iz, giz
  
      write( cnum, fmt="(i3.3)" ) inum
      write( olog, * ) " # Diagnostics of mom."
  
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
        
        if ( ranks == 0 .AND. vel_rank == 0 ) then
  
          write( crank, fmt="(i6.6)" ) ir
          write( srank, fmt="(i1.1)" ) ranks
          open( unit=100000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".phi."//cnum,  &
                status="old", action="read", form="unformatted" )
        !  write( olog, * ) " # Opened file and the unit-ID are ",  &
        !                   "./phi/gkvp_f0.30."//crank//"."//srank//".phi."//cnum, 100000+ir
          open( unit=200000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".Al."//cnum,  &
                status="old", action="read", form="unformatted" )
        !  write( olog, * ) " # Opened file and the unit-ID are ",  &
        !                   "./phi/gkvp_f0.30."//crank//"."//srank//".Al."//cnum, 200000+ir
  
        end if

        if ( vel_rank == 0 ) then

          write( crank, fmt="(i6.6)" ) ir
          write( srank, fmt="(i1.1)" ) ranks
          open( unit=300000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".mom."//cnum,  &
                status="old", action="read", form="unformatted" )
        !  write( olog, * ) " # Opened file and the unit-ID are ",  &
        !                   "./phi/gkvp_f0.30."//crank//"."//srank//".mom."//cnum, 300000+ir

        end if

        if ( zsp_rank == 0 .AND. vel_rank == 0 ) then

          write( crank, fmt="(i6.6)" ) ir
          write( srank, fmt="(i1.1)" ) ranks
          open( unit=400000+ir, file="./phi/gkvp_f0.30."//crank//"."//srank//".trn."//cnum,  &
                status="old", action="read", form="unformatted" )
        !  write( olog, * ) " # Opened file and the unit-ID are ",  &
        !                   "./phi/gkvp_f0.30."//crank//"."//srank//".trn."//cnum, 400000+ir
        end if
  
      end do
  
  
  ! --- time-step loop ---
      loop = 0
      do
  
    if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then            !%%% to restart
  
        if ( flag_phiinkxky .or. flag_phiinxy .or. flag_phiinz  &
             .or. flag_Alinkxky .or. flag_Alinxy .or. flag_Alinz  &
             .or. flag_fluxinkxky .or. flag_fluxinky ) then
           
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
            
            if ( ranks == 0 .AND. vel_rank == 0 ) then
    
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
  
        end if

        if ( flag_fluxinkxky .or. flag_fluxinky ) then

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
  
            if ( vel_rank == 0 ) then
    
              read( unit=300000+ir, iostat=ios ) time, dens, upara, ppara, pperp, qlpara, qlperp
              if ( ios /= 0 ) exit
    
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

        end if

        if ( flag_transinkxky ) then
  
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
  
            if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
    
              read( unit=400000+ir, iostat=ios ) time, entrpy, neint, nmint, dcd
              if ( ios /= 0 ) exit
    
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
  
        end if

    else                             !%%% to restart
  
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
          
          if ( ranks == 0 .AND. vel_rank == 0 ) then
  
            read( unit=100000+ir, iostat=ios )
            if ( ios /= 0 ) exit
            read( unit=200000+ir, iostat=ios )
            if ( ios /= 0 ) exit

          end if

          if ( vel_rank == 0 ) then
            read( unit=300000+ir, iostat=ios )
            if ( ios /= 0 ) exit
  
          end if

          if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
  
            read( unit=400000+ir, iostat=ios )
            if ( ios /= 0 ) exit
  
          end if
  
        end do
  
    end if                           !%%% to restart
  
        if ( ios < 0 ) then
          write( olog, * ) " # All files are completed. iostat, unit-ID = ",  &
                           ios, 100000+ir, 200000+ir, 300000+ir, 400000+ir
          write( olog, * ) ""
          exit
        else if ( ios > 0 ) then
          write( olog, * ) " # I/O error is detected. iostat, unit-ID = ",  &
                           ios, 100000+ir, 200000+ir, 300000+ir, 400000+ir
          write( olog, * ) ""
          exit
        end if
  
    if ( loop >= 0 .and. mod(loop,loopskip)==0 ) then            !%%% to restart

        if ( flag_phiinkxky   ) call phiinkxky( inum, loop, time, gphi )
        if ( flag_phiinxy     ) call phiinxy( inum, loop, time, gphi )
        if ( flag_phiinz      ) call phiinz( inum, loop, time, gphi )
        if ( flag_Alinkxky    ) call Alinkxky( inum, loop, time, gAl )
        if ( flag_Alinxy      ) call Alinxy( inum, loop, time, gAl )
        if ( flag_Alinz       ) call Alinz( inum, loop, time, gAl )
 
        if ( flag_fluxinkxky  ) call fluxinkxky( inum, loop, time, gphi, gAl,  &
                                       gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
        if ( flag_fluxinky    ) call fluxinky( inum, loop, time, gphi, gAl,  &
                                       gdens, gupara, gppara, gpperp, gqlpara, gqlperp )

        if ( flag_transinkxky ) call transinkxky( inum, loop, time, gentrpy, gneint, gnmint, gdcd )
      !  write( olog, * ) " # Data output at loop, time = ", loop, time
        if ( mod(loop,100)==0 ) write(*,*) "inum, loop, time =",inum,loop,time

    end if                           !%%% to restart
  
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
        
        if ( ranks == 0 .AND. vel_rank == 0 ) then
  
          close( unit=100000+ir )
          close( unit=200000+ir )
  
        end if
  
        if ( vel_rank == 0 ) then
  
          close( unit=300000+ir )
  
        end if
  
        if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
  
          close( unit=400000+ir )
  
        end if
  
      end do
  
  
  END SUBROUTINE diag_mom
  
  
!--------------------------------------
  SUBROUTINE phiinkxky( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in kx, ky space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
  ! --- local variables
  
    integer :: mx, my, iz
  
  
      iz = 0
      do my = 0 ,global_ny
        do mx = -nx, nx
          write( ophiinkxky, "(99G17.7E3)" ) kx(mx), gky(my), time, abs( gphi(mx,my,iz) )
        end do
        write( ophiinkxky, * )
      end do
  
  
  END SUBROUTINE phiinkxky
  
  
!--------------------------------------
  SUBROUTINE phiinxy( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in x, y space
  
    use diag_header
    use diag_fft, only : fft_backward_xy
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
  ! --- local variables
  
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: phixy
    integer :: mx, my, iz
  
  
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          phixy(mx,my) = 0._DP
        end do
      end do
  
      iz = 0
      call fft_backward_xy( gphi(:,:,iz), phixy )
  
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( ophiinxy, "(99G17.7E3)" ) xx(mx), yy(my), time, phixy(mx,my)
        end do
        write( ophiinxy, * )
      end do
 
  
  END SUBROUTINE phiinxy
  
  
!--------------------------------------
  SUBROUTINE phiinz( inum, loop, time, gphi )
!--------------------------------------
!     Write phi in zz space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gphi
  
  ! --- local variables
  
    integer :: iconnect, connect_min, connect_max, mxw
    integer :: mx, my, iz
  
  
      mx = 0
      my = 6
  
      if ( dj(my) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( ophiinz, "(99G17.7E3)" ) gzz(iz), time, gphi(mx,my,iz)
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(my) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( ophiinz, "(99G17.7E3)" ) - twopi * real(iconnect) + gzz(iz), time,  &
                                              ck(my)**iconnect * gphi(mxw,my,iz)
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(my) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( ophiinz, "(99G17.7E3)" ) + twopi * real(iconnect) + gzz(iz), time,  &
                                              conjg( ck(my)**iconnect ) * gphi(mxw,my,iz)
            end do
          end do
  
      end if

      write( ophiinz, * )
  
  
  END SUBROUTINE phiinz
  
  
!--------------------------------------
  SUBROUTINE Alinkxky( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in kx, ky space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
  ! --- local variables
  
    integer :: mx, my, iz
  
  
      iz = 0
      do my = 0 ,global_ny
        do mx = -nx, nx
          write( oAlinkxky, "(99G17.7E3)" ) kx(mx), gky(my), time, abs( gAl(mx,my,iz) )
        end do
        write( oAlinkxky, * )
      end do
  
  
  END SUBROUTINE Alinkxky
  
  
!--------------------------------------
  SUBROUTINE Alinxy( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in x, y space
  
    use diag_header
    use diag_fft, only : fft_backward_xy
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
  ! --- local variables
  
    real(kind=DP), dimension(0:2*nxw-1,0:2*nyw-1) :: Alxy
    integer :: mx, my, iz
  
  
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          Alxy(mx,my) = 0._DP
        end do
      end do
  
      iz = 0
      call fft_backward_xy( gAl(:,:,iz), Alxy )
  
      do my = 0, 2*nyw-1
        do mx = 0, 2*nxw-1
          write( oAlinxy, "(99G17.7E3)" ) xx(mx), yy(my), time, Alxy(mx,my)
        end do
        write( oAlinxy, * )
      end do
  
  
  END SUBROUTINE Alinxy
  
  
!--------------------------------------
  SUBROUTINE Alinz( inum, loop, time, gAl )
!--------------------------------------
!     Write Al in zz space
  
    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)                    :: inum
    integer,          intent(in)                    :: loop
    real(kind=DP),    intent(in)                    :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: gAl
  
  ! --- local variables
  
    integer :: iconnect, connect_min, connect_max, mxw
    integer :: mx, my, iz
  
  
      mx = 0
      my = 6
  
      if ( dj(my) == 0 ) then
  
          do iz = -global_nz, global_nz-1
            write( oAlinz, "(99G17.7E3)" ) gzz(iz), time, gAl(mx,my,iz)
          end do
  
      else
  
        connect_min = int( ( nx + mx ) / abs( dj(my) ) )
        if ( connect_min .ne. 0 ) then
          do iconnect = connect_min, 1, -1
            mxw = mx+iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( oAlinz, "(99G17.7E3)" ) - twopi * real(iconnect) + gzz(iz), time,  &
                                             ck(my)**iconnect * gAl(mxw,my,iz)
            end do
          end do
        end if
  
        connect_max = int( ( nx - mx ) / abs( dj(my) ) )
          do iconnect = 0, connect_max
            mxw = mx-iconnect*dj(my)
            do iz = -global_nz, global_nz-1
              write( oAlinz, "(99G17.7E3)" ) + twopi * real(iconnect) + gzz(iz), time,  &
                                             conjg( ck(my)**iconnect ) * gAl(mxw,my,iz)
            end do
          end do
  
      end if

      write( oAlinz, * )
  
  
  END SUBROUTINE Alinz


!--------------------------------------
  SUBROUTINE fluxinkxky( inum, loop, time, gphi, gAl,  &
                         gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!--------------------------------------
!     Write particle and heat fluxes in kx, ky

    use diag_header
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)  :: inum
    integer,          intent(in)  :: loop
    real(kind=DP),    intent(in)  :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1)      :: gphi, gAl
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara,  &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
  
  ! --- local variables
  
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,0:ns-1) :: p_flux_es, p_flux_em, h_flux_es, h_flux_em
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1) :: wc3
    integer :: mx, my, iz, is
  
  
  !--- Electrostatic part of particle fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wc3(mx,my,iz) = conjg( - ui * gky(my) * gphi(mx,my,iz) ) &
                                                   * gdens(mx,my,iz,is)
            end do
          end do
        end do
        call thet_ave_z ( wc3, p_flux_es(:,:,is) )
      end do
  
  
  !--- Magnetic part of particle fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wc3(mx,my,iz) = conjg(   ui * gky(my) * gAl(mx,my,iz) ) &
                                                 * gupara(mx,my,iz,is)
            end do
          end do
        end do
        call thet_ave_z ( wc3, p_flux_em(:,:,is) )
      end do
  
  
  !--- Electrostatic part of heat fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wc3(mx,my,iz) = conjg( - ui * gky(my) * gphi(mx,my,iz) ) &
                         * ( gppara(mx,my,iz,is) + gpperp(mx,my,iz,is) )
            end do
          end do
        end do
        call thet_ave_z ( wc3, h_flux_es(:,:,is) )
      end do
  
  
  !--- Magnetic part of heat fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          do my = 0, global_ny
            do mx = -nx, nx
              wc3(mx,my,iz) = conjg(   ui * gky(my) * gAl(mx,my,iz) ) &
                      * ( gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is) )
            end do
          end do
        end do
        call thet_ave_z ( wc3, h_flux_em(:,:,is) )
      end do
  
      do my = 0, global_ny
        do mx = -nx, nx
          write( ofluxinkxky, "(99G17.7E3)" ) kx(mx), gky(my), time,  &
                            real( p_flux_es(mx,my,0:ns-1), kind=DP ), &
                            real( p_flux_em(mx,my,0:ns-1), kind=DP ), &
                            real( h_flux_es(mx,my,0:ns-1), kind=DP ), &
                            real( h_flux_em(mx,my,0:ns-1), kind=DP )
        end do
        write( ofluxinkxky, * )
      end do
  
  
  END SUBROUTINE fluxinkxky


!--------------------------------------
  SUBROUTINE transinkxky( inum, loop, time, gentrpy, gneint, gnmint, gdcd )
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
  
    integer :: mx, my, iz, is
  
  
      do my = 0, global_ny
        do mx = -nx, nx
          write( otransinkxky, "(99G17.7E3)" ) kx(mx), gky(my), time,  &
                                               gentrpy(mx,my,0:ns-1),  &
                                                gneint(mx,my,0:ns-1),  &
                                                gnmint(mx,my,0:ns-1),  &
                                                  gdcd(mx,my,0:ns-1)
        end do
        write( otransinkxky, * )
      end do
  
  
  END SUBROUTINE transinkxky


!--------------------------------------
  SUBROUTINE fluxinky( inum, loop, time, gphi, gAl,  &
                       gdens, gupara, gppara, gpperp, gqlpara, gqlperp )
!--------------------------------------
!     Write particle and heat fluxes in time, ky

    use diag_header
    use diag_fft, only : fft_backward_x
  
    implicit none
  
  ! --- argument
  
    integer,          intent(in)  :: inum
    integer,          intent(in)  :: loop
    real(kind=DP),    intent(in)  :: time
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1)      :: gphi, gAl
    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: gdens, gupara,  &
                                                                     gppara, gpperp, &
                                                                     gqlpara, gqlperp
  
  ! --- local variables
  
    complex(kind=DP),  &
      dimension(0:2*nxw-1,0:global_ny,-global_nz:global_nz-1) :: phixxky, Alxxky
    complex(kind=DP),  &
      dimension(0:2*nxw-1,0:global_ny,-global_nz:global_nz-1,0:ns-1) :: densxxky, uparaxxky, &
                                                                        presxxky, qparaxxky
    real(kind=DP),  &
      dimension(0:global_ny,0:ns-1) :: p_flux_es, p_flux_em, h_flux_es, h_flux_em,  &
                                       p_ampl_es, p_ampl_em, h_ampl_es, h_ampl_em
    complex(kind=DP),  &
      dimension(-nx:nx,0:global_ny) :: wc2
    real(kind=DP) :: fct
    integer :: mx, my, iz, is
  
  
      p_flux_es(:,:) = 0._DP
      p_flux_em(:,:) = 0._DP
      h_flux_es(:,:) = 0._DP
      h_flux_em(:,:) = 0._DP
      p_ampl_es(:,:) = 0._DP
      p_ampl_em(:,:) = 0._DP
      h_ampl_es(:,:) = 0._DP
      h_ampl_em(:,:) = 0._DP
  
      do iz = -global_nz, global_nz-1
        call fft_backward_x ( gphi(:,:,iz), phixxky(:,:,iz) )
        call fft_backward_x (  gAl(:,:,iz),  Alxxky(:,:,iz) )
      end do
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          call fft_backward_x (  gdens(:,:,iz,is),  densxxky(:,:,iz,is) )
          call fft_backward_x ( gupara(:,:,iz,is), uparaxxky(:,:,iz,is) )
          do my = 0, global_ny
            do mx = -nx, nx
              wc2(mx,my) = gppara(mx,my,iz,is) + gpperp(mx,my,iz,is)
            end do
          end do
          call fft_backward_x ( wc2(:,:),  presxxky(:,:,iz,is) )
          do my = 0, global_ny
            do mx = -nx, nx
              wc2(mx,my) = gqlpara(mx,my,iz,is) + gqlperp(mx,my,iz,is)
            end do
          end do
          call fft_backward_x ( wc2(:,:), qparaxxky(:,:,iz,is) )
        end do
      end do

  !--- Electrostatic part of particle fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / cfsrf / real( 2*nxw, kind=DP )
          do my = 0, global_ny
            do mx = 0, 2*nxw-1
              p_flux_es(my,is) = p_flux_es(my,is) + 2._DP * fct * real(  &
                                 conjg( - ui * gky(my) * phixxky(mx,my,iz) )  &
                                            * densxxky(mx,my,iz,is), kind=DP )
              p_ampl_es(my,is) = p_ampl_es(my,is) + 2._DP * fct * abs(  &
                                 phixxky(mx,my,iz)* densxxky(mx,my,iz,is) )
            end do
          end do
        end do
      end do
  
  
  !--- Magnetic part of particle fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / cfsrf / real( 2*nxw, kind=DP )
          do my = 0, global_ny
            do mx = 0, 2*nxw-1
              p_flux_em(my,is) = p_flux_em(my,is) + 2._DP * fct * real(  &
                                 conjg(   ui * gky(my) * Alxxky(mx,my,iz) )  &
                                          * uparaxxky(mx,my,iz,is), kind=DP )
              p_ampl_em(my,is) = p_ampl_em(my,is) + 2._DP * fct * abs(  &
                                 Alxxky(mx,my,iz)* uparaxxky(mx,my,iz,is) )
            end do
          end do
        end do
      end do
  
  
  !--- Electrostatic part of heat fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / cfsrf / real( 2*nxw, kind=DP )
          do my = 0, global_ny
            do mx = 0, 2*nxw-1
              h_flux_es(my,is) = h_flux_es(my,is) + 2._DP * fct * real(  &
                                 conjg( - ui * gky(my) * phixxky(mx,my,iz) )  &
                                            * presxxky(mx,my,iz,is), kind=DP )
              h_ampl_es(my,is) = h_ampl_es(my,is) + 2._DP * fct * abs(  &
                                 phixxky(mx,my,iz)* presxxky(mx,my,iz,is) )
            end do
          end do
        end do
      end do
  
  
  !--- Magnetic part of heat fluxes ---
      do is = 0, ns-1
        do iz = -global_nz, global_nz-1
          fct = grootg(iz) / cfsrf / real( 2*nxw, kind=DP )
          do my = 0, global_ny
            do mx = 0, 2*nxw-1
              h_flux_em(my,is) = h_flux_em(my,is) + 2._DP * fct * real(  &
                                 conjg(   ui * gky(my) * Alxxky(mx,my,iz) )  &
                                          * qparaxxky(mx,my,iz,is), kind=DP )
              h_ampl_em(my,is) = h_ampl_em(my,is) + 2._DP * fct * abs(  &
                                 Alxxky(mx,my,iz)* qparaxxky(mx,my,iz,is) )
            end do
          end do
        end do
      end do
  
      do my = 0, global_ny
        write( ofluxinky, "(99G17.7E3)" ) gky(my), time,  &
                                   p_flux_es(my,0:ns-1),  &
                                   p_flux_em(my,0:ns-1),  &
                                   h_flux_es(my,0:ns-1),  &
                                   h_flux_em(my,0:ns-1),  &
                                   p_ampl_es(my,0:ns-1),  &
                                   p_ampl_em(my,0:ns-1),  &
                                   h_ampl_es(my,0:ns-1),  &
                                   h_ampl_em(my,0:ns-1)
      end do
      write( ofluxinky, * )
  
  
  END SUBROUTINE fluxinky


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
  
      do iz = -global_nz, global_nz-1
        fct = grootg(iz) / cfsrf
        do my = 0, global_ny
          do mx = -nx, nx
            wa(mx,my)   = wa(mx,my) + fct * wn(mx,my,iz)
          end do
        end do
      end do
  
  
  END SUBROUTINE thet_ave_z


END PROGRAM diag
