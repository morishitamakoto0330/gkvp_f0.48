MODULE reduce_header
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
  logical, parameter :: flag_phi = .true.
!%%%%%%%%%%%%%%%%%%%%%%

!%%% GKV parameters %%%
  integer, parameter :: nxw = 512, nyw = 512
  integer, parameter :: nx = 320, global_ny = 319  ! 2/3 de-aliasing rule
  integer, parameter :: global_nz = 32, global_nv = 48, global_nm = 15
  integer, parameter :: nzb = 2
  integer, parameter :: nprocw = 32, nprocz = 8, nprocv = 12, nprocm = 2, nprocs = 2
  real(kind=DP), parameter ::  vmax = 4._DP
!%%%%%%%%%%%%%%%%%%%%%%

!%%% Data reduction %%%
  integer, parameter :: loopskip = 100
  integer, parameter :: nzskip = 4
  integer, parameter :: rglobal_nz = global_nz/nzskip
  integer, parameter :: rnz = global_nz/nprocz/nzskip
!%%%%%%%%%%%%%%%%%%%%%%

  integer, parameter :: nxw_size = (2*nxw-1)/nprocw
  integer, parameter :: ny       = global_ny / nprocw
  integer, parameter :: nz = global_nz / nprocz,          &
                        nv = global_nv / nprocv,          &
                        nm = (global_nm + 1) / nprocm - 1,&
                        ns = nprocs
  integer, parameter :: nproc = nprocw * nprocz * nprocv * nprocm * nprocs

  real(kind=DP),    parameter :: pi  = 3.141592653589793_DP, twopi = pi * 2._DP
  real(kind=DP),    parameter :: eps = 0.0000000001_DP
  complex(kind=DP), parameter :: ui  = ( 0._DP, 1._DP )

! --- unit numbers for I/O
  integer, parameter :: olog = 10


END MODULE reduce_header


PROGRAM reduce
!-------------------------------------------------------------------------------
!
!     Data diagnostics from binary output
!
!                                   by S.Maeyama  (Dec. 2012)
!
!-------------------------------------------------------------------------------

  use reduce_header

  implicit none

  integer :: inum, loop
  character*3 :: cnum


    loop = 0

    do inum = snum, enum

      if ( flag_phi ) call reduce_phi( inum, loop )

    end do


CONTAINS


!--------------------------------------
  SUBROUTINE reduce_phi( inum, loop )
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
    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: phi
    integer :: rankg, ranks, rank, rankw, rankz, rankv, rankm,  &
               scolor, wcolor, zcolor, vcolor,                  &
               spc_rank, fft_rank, zsp_rank, vel_rank
    integer :: mx, my, iz
  
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
        
        if ( ranks == 0 .and. vel_rank == 0 ) then
  
          write( crank, fmt="(i6.6)" ) ir
          write( srank, fmt="(i1.1)" ) ranks
          open( unit=10, &
                file="../phi/gkvp_f0.30."//crank//"."//srank//".phi."//cnum,  &
                status="old", action="read", form="unformatted" )
          open( unit=20, &
                file="./reduced_phi/gkvp_f0.30."//crank//"."//srank//".phi."//cnum,  &
                action="write", form="unformatted" )

          do

            if ( mod(loop,loopskip) == 0 ) then

              read( unit=10, iostat=ios ) time, phi
              if ( ios /= 0 ) exit

              rphi = phi

              if ( mod(loop,100)==0 ) then
                write(*,*) "inum, loop, time =", inum, loop, time
              end if

            else

              read( unit=10, iostat=ios )
              if ( ios /= 0 ) exit

            end if

            loop = loop + 1

          end do

          close( unit=10 )
          close( unit=20 )

        end if

      end do
  
  
  END SUBROUTINE reduce_phi
  
  
