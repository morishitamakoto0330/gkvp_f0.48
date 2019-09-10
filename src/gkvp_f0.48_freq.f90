MODULE GKV_freq
!-------------------------------------------------------------------------------
!
!     Module for a evaluation of linear growth rate and real frequency
!
!     "freq_set" - set parameters.
!     "freq_reset" - reset parameters.
!     "freq_convcheck" - evaluate omega.
!
!                                              by S. Maeyama (June 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  real(kind=DP) :: time0
  complex(kind=DP), allocatable, save, dimension(:,:) :: phi0
  complex(kind=DP), allocatable, save, dimension(:) :: omega0
  real(kind=DP), allocatable, save, dimension(:) :: phi0_norm2

  complex(kind=DP), allocatable, save, dimension(:) :: omega_g
  complex(kind=DP), allocatable, save, dimension(:) :: diff_g
  real(kind=DP), allocatable, save, dimension(:) :: ineq_g
  logical, save, dimension(0:(ny+1)*nprocw-1) :: freq_conv = .false.

  real(kind=DP), parameter :: eps_omega = 1.d-4, &
                              eps_gamma = 1.d-4, &
                              eps_ineq  = 1.d-6

  public :: freq_write_frq, freq_write_dsp, freq_conv


CONTAINS


!--------------------------------------
  SUBROUTINE freq_set ( time, phi )
!--------------------------------------

    real(kind=DP), intent(in)           :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)   :: phi

    real(kind=DP), dimension(0:ny) :: wr3
    integer :: my, iz


      time0 = time
  
      allocate( phi0(0:ny,-nz:nz-1) )
      allocate( omega0(0:ny) )
      allocate( phi0_norm2(0:ny) )
      allocate( omega_g(0:(ny+1)*nprocw-1) )
      allocate( diff_g(0:(ny+1)*nprocw-1) )
      allocate( ineq_g(0:(ny+1)*nprocw-1) )
      !allocate( freq_conv(0:(ny+1)*nprocw-1) )
  
      phi0(:,:) = ( 0._DP, 0._DP )
      omega0(:) = ( 0._DP, 0._DP )
      phi0_norm2(:) = 0._DP
      omega_g(:) = ( 0._DP, 0._DP )
      diff_g(:) = ( 0._DP, 0._DP )
      ineq_g(:) = 0._DP
      !freq_conv(:) = .false.
      
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          phi0(my,iz) = phi(0,my,iz)
        end do
      end do

      phi0_norm2(:) = 0._DP
      wr3(:) = 0._DP
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          wr3(my) = wr3(my) + abs( phi(0,my,iz) )**2
        end do
      end do
      call MPI_Allreduce( wr3, phi0_norm2, ny+1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_set at t=",time
                                      !%%%%%%%%%%%%%%%%%

  END SUBROUTINE freq_set


!--------------------------------------
  SUBROUTINE freq_reset
!--------------------------------------


      deallocate( phi0 )
      deallocate( omega0 )
      deallocate( phi0_norm2 )
      deallocate( omega_g )
      deallocate( diff_g )
      deallocate( ineq_g )
      !deallocate( freq_conv )
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_reset"
                                      !%%%%%%%%%%%%%%%%%

  END SUBROUTINE freq_reset


!--------------------------------------
  SUBROUTINE freq_write_frq ( time, phi )
!--------------------------------------

    real(kind=DP), intent(in)           :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)   :: phi

    complex(kind=DP), dimension(0:ny) :: phi0phi, omega_l, diff_l
    real(kind=DP), dimension(0:ny) :: phi_norm2, ineq_l
    complex(kind=DP), dimension(0:ny) :: wc3
    real(kind=DP), dimension(0:ny) :: wr3
    integer :: my, iz

    integer, save ::  iflg
    data iflg / 0 /


      if( iflg == 0 ) then
        iflg = 1
        call freq_set ( time, phi )
        return
      end if


!- calculate interior products -
      phi0phi(:) = (0._DP, 0._DP)
      wc3(:) = (0._DP, 0._DP)
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          wc3(my) = wc3(my) + conjg( phi0(my,iz) ) * phi(0,my,iz)
        end do
      end do
      call MPI_Allreduce( wc3, phi0phi, ny+1, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

      phi_norm2(:) = 0._DP
      wr3(:) = 0._DP
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          wr3(my) = wr3(my) + abs( phi(0,my,iz) )**2
        end do
      end do
      call MPI_Allreduce( wr3, phi_norm2, ny+1, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

!- calculate frequency -
      do my = ist1_y, iend_y
        omega_l(my) = log( phi0phi(my) / phi0_norm2(my) )  &
                                    / ( ui * ( time0 - time ) )
      end do
      call MPI_Allgather( omega_l, ny+1, MPI_DOUBLE_COMPLEX, &
                          omega_g, ny+1, MPI_DOUBLE_COMPLEX, &
                          fft_comm_world, ierr_mpi )

!- convergence check -
      do my = ist1_y, iend_y
        diff_l(my) = abs(real(omega_l(my) - omega0(my), kind=DP)  &
                                    / real(omega_l(my), kind=DP)) &
                   + ui * abs(aimag(omega_l(my) - omega0(my))     &
                                    / aimag(omega_l(my)) )
      end do
      call MPI_Allgather( diff_l, ny+1, MPI_DOUBLE_COMPLEX, &
                          diff_g, ny+1, MPI_DOUBLE_COMPLEX, &
                          fft_comm_world, ierr_mpi )

      do my = ist1_y, iend_y
        ineq_l(my) = abs(phi0phi(my))**2 / (phi0_norm2(my) * phi_norm2(my))
      end do
      call MPI_Allgather( ineq_l, ny+1, MPI_DOUBLE_PRECISION, &
                          ineq_g, ny+1, MPI_DOUBLE_PRECISION, &
                          fft_comm_world, ierr_mpi )

      do my = 1, global_ny
        if ( real( diff_g(my), kind=DP ) < eps_omega .and.  &
             aimag( diff_g(my) ) < eps_gamma .and.          &
             (1._DP - ineq_g(my)) < eps_ineq ) then
          freq_conv(my) = .true.
        else
          freq_conv(my) = .false.
        end if
      end do

      freq_conv(0) = .true.
      if ( global_ny < (ny+1)*nprocw-1 ) then
        do my = global_ny+1, (ny+1)*nprocw-1
          diff_g(my) = ( 0._DP, 0._DP )
          omega_g(my) = ( 0._DP, 0._DP )
          ineq_g(my) = 0._DP
          freq_conv(my) = .true.
        end do
      end if

!- remember the values -
      time0 = time
      do iz = -nz, nz-1
        do my = ist1_y, iend_y
          phi0(my,iz) = phi(0,my,iz)
        end do
      end do
      do my = ist1_y, iend_y
        omega0(my) = omega_l(my)
        phi0_norm2(my) = phi_norm2(my)
      end do

!- write hst/*.frq.* -
      if ( rankg == 0 ) then
        write( ofrq, '(9999G17.7e3)' ) time, (omega_g(my), my=1,global_ny)
      end if
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_write_frq at t=",time
                                      !%%%%%%%%%%%%%%%%%

END SUBROUTINE freq_write_frq


!--------------------------------------
  SUBROUTINE freq_write_dsp
!--------------------------------------

    integer :: my


      if ( rankg == 0 ) then
        write( odsp, '(99A17)' ) "#              ky","frequency","growthrate",&
                                           "diff(freq)","diff(grow)","1-ineq"
        do my = 1, global_ny
          if ( freq_conv(my) ) then
            write( odsp, '(9999G17.7e3)' ) ky(1) * real( my, kind=DP ),       &
                       real( omega_g(my), kind=DP ), aimag( omega_g(my) ),    &
                       real( diff_g(my), kind=DP ), aimag( diff_g(my) ),      &
                       1._DP - ineq_g(my)
          else
            write( odsp, '(A2,9999G17.7e3)' ) "# ", ky(1) * real( my, kind=DP ),&
                       real( omega_g(my), kind=DP ), aimag( omega_g(my) ),      &
                       real( diff_g(my), kind=DP ), aimag( diff_g(my) ),        &
                       1._DP - ineq_g(my)
          end if
        end do
      end if
  
      call freq_reset
                                      !%%% For debug %%%
                                      !write(olog,*) "freq_write_dsp"
                                      !%%%%%%%%%%%%%%%%%


  END SUBROUTINE freq_write_dsp


END MODULE GKV_freq
