MODULE GKV_trans
!-------------------------------------------------------------------------------
!
!    Entropy transfer diagnostics
!
!      GKV-plus f0.26 ( S.Maeyama, Nov 2012 )
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fft,   only: fft_pre,  &
           fft_backward_Xfft, fft_backward_chXY, fft_backward_Yfft, &
           fft_forward_Yfft, fft_forward_chYX, fft_forward_Xfft
  use GKV_intgrl, only: intgrl_thet, intgrl_v0_moment
  use GKV_clock, only: clock_sta, clock_end

  implicit none

  private

  ! variables divided in Y
  complex(kind=DP), save, &
    dimension(0:2*nxw-1,0:ny,-nz:nz-1,1:2*nv,0:nm) :: exbdf1, exbdf2
  complex(kind=DP), save, &
    dimension(0:2*nxw-1,0:ny,-nz:nz-1,0:nm)        :: exw, eyw, bxw, byw
  complex(kind=DP), save, &
    dimension(0:2*nxw-1,0:ny)                 :: uikx, uiky

  ! variables divided in X
  complex(kind=DP), save, &
    dimension(0:nyw,0:nxw_size,-nz:nz-1,1:2*nv,0:nm) :: exbdf1_xw, exbdf2_xw   ! X,Y order changed
  complex(kind=DP), save, &
    dimension(0:nyw,0:nxw_size,-nz:nz-1,0:nm)        :: exw_xw, eyw_xw, bxw_xw, byw_xw   ! X,Y order changed

  integer, save, &
    dimension(:), allocatable :: triad_diag_mxt, triad_diag_myt
  integer, save :: nbuff

  public   trans_sum, trans_triad


CONTAINS


!--------------------------------------
  SUBROUTINE trans_sum ( ff, phi, Al, neint, nmint )
!--------------------------------------
!     Check the entropy balance equation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)                :: neint, nmint

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: wf, wef, wbf
    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi
    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: wc3
    complex(kind=DP), dimension(-nx:nx,0:ny) :: wc2
    integer  ::  mx, my, iz, iv, im


      allocate( wf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( wef(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( wbf(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( psi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )
      allocate( chi(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) )

!$OMP parallel do
      do im = 0, nm
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                psi(mx,my,iz,im)   = j0(mx,my,iz,im) * phi(mx,my,iz)
                chi(mx,my,iz,im)   = j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
      end do

!$OMP parallel do
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wf(mx,my,iz,iv,im)   = ff(mx,my,iz,iv,im) + sgn(ranks) * Znum(ranks) &
                          * fmx(iz,iv,im) / tau(ranks) * psi(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      call trans_NL_term_em( wf, psi, chi, wef, wbf )

!$OMP parallel do
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wef(mx,my,iz,iv,im)   = - fcs(ranks) / Znum(ranks) * wef(mx,my,iz,iv,im)  &
                      * tau(ranks) * conjg( wf(mx,my,iz,iv,im) ) / fmx(iz,iv,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment ( wef, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          neint(mx,my) = real( wc2(mx,my), kind=DP )
        end do
      end do


!$OMP parallel do
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wbf(mx,my,iz,iv,im)   = - fcs(ranks) / Znum(ranks) * wbf(mx,my,iz,iv,im)  &
                      * tau(ranks) * conjg( wf(mx,my,iz,iv,im) ) / fmx(iz,iv,im)
              end do
            end do
          end do
        end do
      end do

      call intgrl_v0_moment ( wbf, wc3 )

      call intgrl_thet ( wc3, wc2 )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          nmint(mx,my) = real( wc2(mx,my), kind=DP )
        end do
      end do


      deallocate( wf )
      deallocate( wef )
      deallocate( wbf )
      deallocate( psi )
      deallocate( chi )


  END SUBROUTINE trans_sum


!--------------------------------------
  SUBROUTINE trans_NL_term_em( hh, psi, chi, ef, bf )
!--------------------------------------
!  ExB nonlinear term calculation 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)                :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)                 :: ef, bf

    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz)*(2*nv),0:nprocw-1,0:nm) :: send_df1, recv_df1
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz)*(2*nv),0:nprocw-1,0:nm) :: send_df2, recv_df2
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_exw, recv_exw
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_eyw, recv_eyw
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_bxw, recv_bxw
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_byw, recv_byw
    integer       ::  num_trans_exbdf, num_trans_exyw, num_trans_bxyw
    integer, save ::  iflg
    integer       ::  mx, my, im

    data iflg / 0 /

                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_input",1410,1)
      if( iflg == 0 ) then

        iflg = 1

        uikx(:,:) = ( 0._DP, 0._DP )
        uiky(:,:) = ( 0._DP, 0._DP )

!$OMP parallel do
        do my = ist_y, iend_y
          do mx = 0, nx
            uikx(mx,my) = kx(mx) * ui
            uiky(mx,my) = ky(my) * ui
          end do
        end do

!$OMP parallel do
        do my = ist_y, iend_y
          do mx = -nx, -1
            uikx(2*nxw+mx,my) = kx(mx) * ui
            uiky(2*nxw+mx,my) = ky(my) * ui
          end do
        end do

      end if
                                         ! call fapp_stop("nlterm_input",1410,1)
                                           call clock_end(1410)

      num_trans_exbdf = (2*nz)*(2*nv)
      num_trans_exyw  = (2*nz)
      num_trans_bxyw  = (2*nz)

!$OMP parallel default (none)                     &
!$OMP shared (hh,psi,chi,ef,bf)                   &
!$OMP shared (exbdf1,send_df1,recv_df1,exbdf1_xw) &
!$OMP shared (exbdf2,send_df2,recv_df2,exbdf2_xw) &
!$OMP shared (   exw,send_exw,recv_exw,   exw_xw) &
!$OMP shared (   eyw,send_eyw,recv_eyw,   eyw_xw) &
!$OMP shared (   bxw,send_bxw,recv_bxw,   bxw_xw) &
!$OMP shared (   byw,send_byw,recv_byw,   byw_xw) &
!$OMP shared (num_trans_exbdf,num_trans_exyw,num_trans_bxyw) &
!$OMP private (im)

! -------------------------------------------------------------------
!     Backward FFT region
!     5 routines are overlaped:
!       trans_sum_input, backward FFT in X, MPI_AlltoAll from X to Y,  
!       backward FFT in Y, real space cal.
! -------------------------------------------------------------------
      call trans_sum_input (     hh(:,:,:,:,0), psi(:,:,:,0), chi(:,:,:,0),  &
                       exbdf1(:,:,:,:,0), exbdf2(:,:,:,:,0),           &
                            exw(:,:,:,0),      eyw(:,:,:,0),           &
                            bxw(:,:,:,0),      byw(:,:,:,0)           )
!$OMP barrier

      call trans_sum_input (     hh(:,:,:,:,1), psi(:,:,:,1), chi(:,:,:,1),  &
                       exbdf1(:,:,:,:,1), exbdf2(:,:,:,:,1),           &
                            exw(:,:,:,1),      eyw(:,:,:,1),           &
                            bxw(:,:,:,1),      byw(:,:,:,1)           )

      call fft_backward_Xfft ( exbdf1(:,:,:,:,0), send_df1(:,:,:,:,0), num_trans_exbdf )
      call fft_backward_Xfft ( exbdf2(:,:,:,:,0), send_df2(:,:,:,:,0), num_trans_exbdf )
      call fft_backward_Xfft (      exw(:,:,:,0), send_exw(:,:,:,:,0), num_trans_exyw  )
      call fft_backward_Xfft (      eyw(:,:,:,0), send_eyw(:,:,:,:,0), num_trans_exyw  )
      call fft_backward_Xfft (      bxw(:,:,:,0), send_bxw(:,:,:,:,0), num_trans_bxyw  )
      call fft_backward_Xfft (      byw(:,:,:,0), send_byw(:,:,:,:,0), num_trans_bxyw  )
!$OMP barrier

      im = 0
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_bxw(:,:,:,:,im), recv_bxw(:,:,:,:,im), num_trans_bxyw  )
        call fft_backward_chXY ( send_byw(:,:,:,:,im), recv_byw(:,:,:,:,im), num_trans_bxyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call trans_sum_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      bxw(:,:,:,im+1), send_bxw(:,:,:,:,im+1), num_trans_bxyw  )
        call fft_backward_Xfft (      byw(:,:,:,im+1), send_byw(:,:,:,:,im+1), num_trans_bxyw  )
!$OMP barrier

      im = 1
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_bxw(:,:,:,:,im), recv_bxw(:,:,:,:,im), num_trans_bxyw  )
        call fft_backward_chXY ( send_byw(:,:,:,:,im), recv_byw(:,:,:,:,im), num_trans_bxyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call trans_sum_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      bxw(:,:,:,im+1), send_bxw(:,:,:,:,im+1), num_trans_bxyw  )
        call fft_backward_Xfft (      byw(:,:,:,im+1), send_byw(:,:,:,:,im+1), num_trans_bxyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,im-1),      bxw_xw(:,:,:,im-1), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,im-1),      byw_xw(:,:,:,im-1), num_trans_bxyw  )
!$OMP barrier

      do im = 2, nm-2
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_bxw(:,:,:,:,im), recv_bxw(:,:,:,:,im), num_trans_bxyw  )
        call fft_backward_chXY ( send_byw(:,:,:,:,im), recv_byw(:,:,:,:,im), num_trans_bxyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call trans_sum_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      bxw(:,:,:,im+1), send_bxw(:,:,:,:,im+1), num_trans_bxyw  )
        call fft_backward_Xfft (      byw(:,:,:,im+1), send_byw(:,:,:,:,im+1), num_trans_bxyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,im-1),      bxw_xw(:,:,:,im-1), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,im-1),      byw_xw(:,:,:,im-1), num_trans_bxyw  )

        call trans_sum_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
                                  exw_xw(:,:,:,im-2),      eyw_xw(:,:,:,im-2), & 
                                  bxw_xw(:,:,:,im-2),      byw_xw(:,:,:,im-2) )
!$OMP barrier
      end do

      im = nm-1
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_bxw(:,:,:,:,im), recv_bxw(:,:,:,:,im), num_trans_bxyw  )
        call fft_backward_chXY ( send_byw(:,:,:,:,im), recv_byw(:,:,:,:,im), num_trans_bxyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      bxw(:,:,:,im+1), send_bxw(:,:,:,:,im+1), num_trans_bxyw  )
        call fft_backward_Xfft (      byw(:,:,:,im+1), send_byw(:,:,:,:,im+1), num_trans_bxyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,im-1),      bxw_xw(:,:,:,im-1), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,im-1),      byw_xw(:,:,:,im-1), num_trans_bxyw  )

        call trans_sum_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
                                  exw_xw(:,:,:,im-2),      eyw_xw(:,:,:,im-2), & 
                                  bxw_xw(:,:,:,im-2),      byw_xw(:,:,:,im-2) )
!$OMP barrier


! -------------------------------------------------------------------
!     Backward & forward FFT region
!     6 routines are overlaped: 
!       MPI_AlltoAll from X to Y, backward FFT in Y, real space cal.,
!       forward FFT in Y, MPI_AlltoAll form Y to X, forward FFT in X
! -------------------------------------------------------------------
      im = nm
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_bxw(:,:,:,:,im), recv_bxw(:,:,:,:,im), num_trans_bxyw  )
        call fft_backward_chXY ( send_byw(:,:,:,:,im), recv_byw(:,:,:,:,im), num_trans_bxyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,im-1),      bxw_xw(:,:,:,im-1), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,im-1),      byw_xw(:,:,:,im-1), num_trans_bxyw  )

        call trans_sum_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
                                  exw_xw(:,:,:,im-2),      eyw_xw(:,:,:,im-2), & 
                                  bxw_xw(:,:,:,im-2),      byw_xw(:,:,:,im-2) )

        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,0), send_df1(:,:,:,:,0), num_trans_exbdf )
        call fft_forward_Yfft ( exbdf2_xw(:,:,:,:,0), send_df2(:,:,:,:,0), num_trans_exbdf )
!$OMP barrier

      im = 0
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_forward_chYX ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_forward_Yfft ( exbdf2_xw(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )

        call fft_backward_Yfft ( recv_df1(:,:,:,:,nm), exbdf1_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,nm),      exw_xw(:,:,:,nm), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,nm),      eyw_xw(:,:,:,nm), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,nm),      bxw_xw(:,:,:,nm), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,nm),      byw_xw(:,:,:,nm), num_trans_bxyw  )

        call trans_sum_realspcal ( exbdf1_xw(:,:,:,:,nm-1), exbdf2_xw(:,:,:,:,nm-1), &
                                  exw_xw(:,:,:,nm-1),      eyw_xw(:,:,:,nm-1), & 
                                  bxw_xw(:,:,:,nm-1),      byw_xw(:,:,:,nm-1) )
!$OMP barrier

      im = 1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_forward_chYX ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_forward_Yfft ( exbdf2_xw(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )
        call fft_forward_Xfft ( recv_df2(:,:,:,:,im-1), exbdf2(:,:,:,:,im-1), num_trans_exbdf )

        call trans_sum_realspcal ( exbdf1_xw(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), &
                                  exw_xw(:,:,:,nm),      eyw_xw(:,:,:,nm), &
                                  bxw_xw(:,:,:,nm),      byw_xw(:,:,:,nm) )
!$OMP barrier


! -------------------------------------------------------------------
!     Forward FFT region
!     4 routines are overlaped:
!       forward FFT in Y, MPI_AlltoAll form Y to X,
!       forward FFT in X, trans_sum_output
! -------------------------------------------------------------------
      do im = 2, nm-1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_forward_chYX ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_forward_Yfft ( exbdf2_xw(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )
        call fft_forward_Xfft ( recv_df2(:,:,:,:,im-1), exbdf2(:,:,:,:,im-1), num_trans_exbdf )

        call trans_sum_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
        call trans_sum_output ( exbdf2(:,:,:,:,im-2), bf(:,:,:,:,im-2) )
!$OMP barrier
      end do

      im = nm
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_forward_chYX ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )
        call fft_forward_Xfft ( recv_df2(:,:,:,:,im-1), exbdf2(:,:,:,:,im-1), num_trans_exbdf )

        call trans_sum_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
        call trans_sum_output ( exbdf2(:,:,:,:,im-2), bf(:,:,:,:,im-2) )
!$OMP barrier

      call fft_forward_Xfft ( recv_df1(:,:,:,:,nm), exbdf1(:,:,:,:,nm), num_trans_exbdf )
      call fft_forward_Xfft ( recv_df2(:,:,:,:,nm), exbdf2(:,:,:,:,nm), num_trans_exbdf )

      call trans_sum_output ( exbdf1(:,:,:,:,nm-1), ef(:,:,:,:,nm-1) )
      call trans_sum_output ( exbdf2(:,:,:,:,nm-1), bf(:,:,:,:,nm-1) )
!$OMP barrier

      call trans_sum_output ( exbdf1(:,:,:,:,nm), ef(:,:,:,:,nm) )
      call trans_sum_output ( exbdf2(:,:,:,:,nm), bf(:,:,:,:,nm) )
!$OMP end parallel


  END SUBROUTINE trans_NL_term_em


!--------------------------------------
  SUBROUTINE trans_sum_input ( hh, psi, chi, wkdf1, wkdf2, wkexw, wkeyw, wkbxw, wkbyw )
!--------------------------------------
!     Data input for E x B term calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb)            :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:ny,-nz:nz-1,1:2*nv) :: wkdf1, wkdf2
    complex(kind=DP), intent(out), &
      dimension(0:2*nxw-1,0:ny,-nz:nz-1)        :: wkexw, wkeyw, wkbxw, wkbyw

    complex(kind=DP) :: ww
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1410)
                                         ! call fapp_start("nlterm_input",1410,1)
!$OMP end master

! --- for hh

!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1

          do my = ist_y, iend_y
            do mx = 0, nx
              wkdf1(mx,my,iz,iv) = hh(mx,my,iz,iv)
            end do
          end do

          do my = ist_y, iend_y
            do mx = nx+1, 2*nxw-nx-1
              wkdf1(mx,my,iz,iv) = ( 0._DP, 0._DP )
            end do
          end do

          do my = ist_y, iend_y
            do mx = -nx, -1
              wkdf1(2*nxw+mx,my,iz,iv) = hh(mx,my,iz,iv)
            end do
          end do

          do my = ist_y, iend_y
            do mx = 0, 2*nxw-1
              ww                 = wkdf1(mx,my,iz,iv)
              wkdf1(mx,my,iz,iv) = uikx(mx,my) * ww
              wkdf2(mx,my,iz,iv) = uiky(mx,my) * ww
            end do
          end do

          ! --- reality condition
          !     gradients of (0,0) mode should be zero 
          if( rankw == 0 )  then
            my = 0
              do mx = 1, nx
                wkdf1(2*nxw-mx,my,iz,iv) = conjg( wkdf1(mx,my,iz,iv) )
                wkdf2(2*nxw-mx,my,iz,iv) = conjg( wkdf2(mx,my,iz,iv) )
              end do
          end if

        end do
      end do
!$OMP end do nowait


! --- for ex and ey

!$OMP do schedule (dynamic)
      do iz = -nz, nz-1

        do my = ist_y, iend_y
          do mx = 0, nx
            wkexw(mx,my,iz) = psi(mx,my,iz)
          end do
        end do

        do my = ist_y, iend_y
          do mx = nx+1, 2*nxw-nx-1
            wkexw(mx,my,iz) = ( 0._DP, 0._DP )
          end do
        end do

        do my = ist_y, iend_y
          do mx = -nx, -1
            wkexw(2*nxw+mx,my,iz) = psi(mx,my,iz)
          end do
        end do

        do my = ist_y, iend_y
          do mx = 0, 2*nxw-1
            ww              = wkexw(mx,my,iz)
            wkexw(mx,my,iz) = - uikx(mx,my) * ww
            wkeyw(mx,my,iz) = - uiky(mx,my) * ww
          end do
        end do

        ! --- reality condition
        !     gradients of (0,0) mode should be zero 
        if( rankw == 0 )  then
          my = 0
            do mx = 1, nx
              wkexw(2*nxw-mx,my,iz) = conjg( wkexw(mx,my,iz) )
              wkeyw(2*nxw-mx,my,iz) = conjg( wkeyw(mx,my,iz) )
            end do
        end if

      end do
!$OMP end do nowait


! --- for mx and my

!$OMP do schedule (dynamic)
      do iz = -nz, nz-1

        do my = ist_y, iend_y
          do mx = 0, nx
            wkbxw(mx,my,iz) = chi(mx,my,iz)
          end do
        end do

        do my = ist_y, iend_y
          do mx = nx+1, 2*nxw-nx-1
            wkbxw(mx,my,iz) = ( 0._DP, 0._DP )
          end do
        end do

        do my = ist_y, iend_y
          do mx = -nx, -1
            wkbxw(2*nxw+mx,my,iz) = chi(mx,my,iz)
          end do
        end do

        do my = ist_y, iend_y
          do mx = 0, 2*nxw-1
            ww              = wkbxw(mx,my,iz)
            wkbxw(mx,my,iz) = - uikx(mx,my) * ww
            wkbyw(mx,my,iz) = - uiky(mx,my) * ww
          end do
        end do

        ! --- reality condition
        !     gradients of (0,0) mode should be zero 
        if( rankw == 0 )  then
          my = 0
            do mx = 1, nx
              wkbxw(2*nxw-mx,my,iz) = conjg( wkbxw(mx,my,iz) )
              wkbyw(2*nxw-mx,my,iz) = conjg( wkbyw(mx,my,iz) )
            end do
        end if

      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_input",1410,1)
                                           call clock_end(1410)
!$OMP end master


  END SUBROUTINE trans_sum_input


!--------------------------------------
  SUBROUTINE trans_sum_realspcal ( wkdf1_xw, wkdf2_xw, wkexw_xw, wkeyw_xw, wkbxw_xw, wkbyw_xw )
!--------------------------------------
!     Calculate E x B term in real space

    complex(kind=DP), intent(inout), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1,1:2*nv) :: wkdf1_xw
    complex(kind=DP), intent(inout), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1,1:2*nv) :: wkdf2_xw
    complex(kind=DP), intent(in), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1)        :: wkexw_xw, wkeyw_xw
    complex(kind=DP), intent(in), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1)        :: wkbxw_xw, wkbyw_xw
    complex(kind=DP) :: wk
    real(kind=DP) ::  cef, cs1
    integer  ::  mx, my, iz, iv

      cef = 1._DP / real( nnx*nny, kind=DP )
      cs1 = sqrt( tau(ranks) / Anum(ranks) )

! --- Real space calculation
!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do mx = ist_xw, iend_xw
            do my = 0, nyw
              wk = wkdf1_xw(my,mx,iz,iv)
              wkdf1_xw(my,mx,iz,iv) =                           &
                 cmplx(      real ( wk                   , kind=DP ) &
                           * real ( wkeyw_xw(my,mx,iz),    kind=DP ) &
                       -     real ( wkdf2_xw(my,mx,iz,iv), kind=DP ) &
                           * real ( wkexw_xw(my,mx,iz),    kind=DP ) &
                       ,    aimag ( wk                             ) &
                          * aimag ( wkeyw_xw(my,mx,iz)             ) &
                      -     aimag ( wkdf2_xw(my,mx,iz,iv)          ) &
                          * aimag ( wkexw_xw(my,mx,iz)             ) &
                      , kind=DP ) * cef
              wkdf2_xw(my,mx,iz,iv) =                           &
                 cmplx(      real ( wk                   , kind=DP ) &
                           * real ( - cs1 * vl(iv) * wkbyw_xw(my,mx,iz),    kind=DP ) &
                       -     real ( wkdf2_xw(my,mx,iz,iv), kind=DP ) &
                           * real ( - cs1 * vl(iv) * wkbxw_xw(my,mx,iz),    kind=DP ) &
                       ,    aimag ( wk                             ) &
                          * aimag ( - cs1 * vl(iv) * wkbyw_xw(my,mx,iz)             ) &
                      -     aimag ( wkdf2_xw(my,mx,iz,iv)          ) &
                          * aimag ( - cs1 * vl(iv) * wkbxw_xw(my,mx,iz)             ) &
                      , kind=DP ) * cef
            end do
          end do
        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master


  END SUBROUTINE trans_sum_realspcal


!--------------------------------------
  SUBROUTINE trans_sum_output ( wkdf1, wf )
!--------------------------------------
!     Data output from E x B term calculation

    complex(kind=DP), intent(in), &
      dimension(0:2*nxw-1,0:ny,-nz:nz-1,1:2*nv) :: wkdf1
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)    :: wf

    integer  ::  mx, my, iz, iv

! --- Wavenumber space substitution
!$OMP master
                                           call clock_sta(1450)
                                         ! call fapp_start("nlterm_output",1450,1)
!$OMP end master
!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1

          do my = ist_y, iend_y
            do mx = 0, nx
              wf(mx,my,iz,iv) = wkdf1(mx,my,iz,iv) * baxfactor
            end do
          end do

          do my = ist_y, iend_y
            do mx = -nx, -1
              wf(mx,my,iz,iv) = wkdf1(2*nxw+mx,my,iz,iv) * baxfactor
            end do
          end do

        end do
      end do
!$OMP end do nowait
!$OMP master
                                         ! call fapp_stop("nlterm_output",1450,1)
                                           call clock_end(1450)
!$OMP end master


  END SUBROUTINE trans_sum_output


!--------------------------------------
  SUBROUTINE trans_triad ( time, ff, phi, Al )
!--------------------------------------
!   Triad transfer    T [(kx,ky)|(px,py),(qx,qy)] * delta(kx+px+qx=0,ky+py+qy=0)
!   Setting (kx,ky)=(diagx,diagy) and (qx=-px-kx,qy=-py-ky), 
!   T=T(px,py) represents transfer from (px,py) and (-px-diagx,-py-diagy) to (diagx,diagy).

    real(kind=DP), intent(in)                                     :: time
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                             :: phi, Al

    complex(kind=DP), dimension(:,:,:,:,:), allocatable :: gg
    complex(kind=DP), dimension(:,:,:,:), allocatable :: psi, chi, wg
    complex(kind=DP), dimension(:,:,:), allocatable :: wps, wch
    real(kind=DP), dimension(:,:), allocatable :: jkpq_es, jpqk_es, jqkp_es, &
                                                  jkpq_em, jpqk_em, jqkp_em
    real(kind=DP) :: ceff
    integer  ::  mx, my, iz, iv, im, it, mxt, myt
    character(1)   :: srank
    character(3)   :: cnew
    character(4)   :: cmx, cmy
    character(256)   :: env_string
    integer, save ::  iflg

    data iflg / 0 /
    namelist /triad/ mxt, myt


      !%%% Initialize triad_diag_mxt, triad_diag_myt, nbuff %%%
      if( iflg == 0 ) then

        iflg = 1

        close(inml)
        call getenv ( 'fu05',env_string )
        open(inml, file=env_string )

        allocate(triad_diag_mxt(0:num_triad_diag-1))
        allocate(triad_diag_myt(0:num_triad_diag-1))

        do it = 0, num_triad_diag-1
          read(inml, nml=triad)
          triad_diag_mxt(it) = mxt
          triad_diag_myt(it) = myt

          if ( rank == 0 ) then
            write( srank, fmt="(i1.1)" ) ranks
            write( cnew,  fmt="(i3.3)" ) inum
            write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
            write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
            open( otri, file=trim(f_phi)//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, & 
                        form="unformatted", status="replace" )
            close( otri )
          end if

        end do

        if ( mod(2*nz*(nm+1),nprocw) == 0 ) then
          nbuff = 2*nz*(nm+1)/nprocw
        else
          nbuff = 2*nz*(nm+1)/nprocw + 1
        end if

      end if

      !%%% Set gg, psi, chi %%%
      allocate( gg(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) )
      allocate( psi(-nx:nx,0:ny,-nz:nz-1,0:nm) )
      allocate( chi(-nx:nx,0:ny,-nz:nz-1,0:nm) )
      allocate( wg(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) )
      allocate( wps(-nx:nx,0:global_ny,0:nbuff-1) )
      allocate( wch(-nx:nx,0:global_ny,0:nbuff-1) )
      allocate( jkpq_es(-nx:nx,-global_ny:global_ny) )
      allocate( jpqk_es(-nx:nx,-global_ny:global_ny) )
      allocate( jqkp_es(-nx:nx,-global_ny:global_ny) )
      allocate( jkpq_em(-nx:nx,-global_ny:global_ny) )
      allocate( jpqk_em(-nx:nx,-global_ny:global_ny) )
      allocate( jqkp_em(-nx:nx,-global_ny:global_ny) )
!$OMP parallel workshare
      gg(:,:,:,:,:) = (0._DP, 0._DP)
      psi(:,:,:,:) = (0._DP, 0._DP)
      chi(:,:,:,:) = (0._DP, 0._DP)
!$OMP end parallel workshare

!$OMP parallel do
      do im = 0, nm
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                psi(mx,my,iz,im) = j0(mx,my,iz,im) * phi(mx,my,iz)
                chi(mx,my,iz,im) = j0(mx,my,iz,im) * Al(mx,my,iz)
              end do
            end do
          end do
      end do
!$OMP parallel do
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                gg(mx,my,iz,iv,im) = ff(mx,my,iz,iv,im)                       &
                                   + sgn(ranks) * Znum(ranks) * fmx(iz,iv,im) &
                                              / tau(ranks) * psi(mx,my,iz,im)
              end do
            end do
          end do
        end do
      end do

      !%%% Multiply the jacobian on gg for zz,vl,mu integration %%%
      if ( rankm == 0 ) then
        im = 0
!$OMP parallel do private(mx,my,iz,iv,ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = 0._DP
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        im = 1
!$OMP parallel do private(mx,my,iz,iv,ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi * (     &
                     1._DP + 1._DP/12._DP + 22._DP / 720._DP )&
                   * rootg(iz) / cfsrf                        &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        im = 2
!$OMP parallel do private(mx,my,iz,iv,ceff)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi * ( &
                     1._DP - 11._DP / 720._DP )           &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
!$OMP parallel do private(mx,my,iz,iv,im,ceff)
        do im = 3, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi     &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do
      else
!$OMP parallel do private(mx,my,iz,iv,im,ceff)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              ceff = vp(iz,im) * dvp(iz) * dv * twopi     &
                   * rootg(iz) / cfsrf                    &
                   * tau(ranks) / (2._DP * fmx(iz,iv,im))
              do my = ist_y, iend_y
                do mx = -nx, nx
                  gg(mx,my,iz,iv,im) = sqrt(ceff) * gg(mx,my,iz,iv,im)
                end do
              end do
            end do
          end do
        end do
      end if

      !%%% Transpose from (y,z,v,m) to (zm,z,v,m) decomposition %%%
      call trans_triad_transpose(gg, psi, chi, wg, wps, wch)

      do it = 0, num_triad_diag-1

        !%%% Calculate traid coupling %%%
        call trans_triad_coupling(it, wg, wps, wch, jkpq_es, jpqk_es, jqkp_es,&
                                                    jkpq_em, jpqk_em, jqkp_em)

        !%%% Output %%%
        if ( rank == 0 ) then
          write( srank, fmt="(i1.1)" ) ranks
          write( cnew,  fmt="(i3.3)" ) inum
          write( cmx,   fmt="(i4.4)" ) triad_diag_mxt(it)
          write( cmy,   fmt="(i4.4)" ) triad_diag_myt(it)
          open( otri, file=trim(f_phi)//"s"//srank//"mx"//cmx//"my"//cmy//".tri."//cnew, & 
                      form="unformatted", status="unknown", position="append" )
          write( unit=otri ) time, jkpq_es, jpqk_es, jqkp_es, jkpq_em, jpqk_em, jqkp_em
          close( otri )
        end if

      end do

      deallocate( gg )
      deallocate( psi )
      deallocate( chi )
      deallocate( wg )
      deallocate( wps )
      deallocate( wch )
      deallocate( jkpq_es )
      deallocate( jpqk_es )
      deallocate( jqkp_es )
      deallocate( jkpq_em )
      deallocate( jpqk_em )
      deallocate( jqkp_em )


  END SUBROUTINE trans_triad


!--------------------------------------
  SUBROUTINE trans_triad_transpose ( gg, psi, chi, wg, wps, wch )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: gg
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:nm)        :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) :: wg
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:global_ny,0:nbuff-1)        :: wps, wch

    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,1:2*nv,0:nbuff-1,0:nprocw-1) :: send_gg, recv_gg
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,0:nbuff-1,0:nprocw-1) :: send_psi, recv_psi, &
                                                     send_chi, recv_chi
    integer  ::  mx, my, iz, iv, im, izm, ibuff, iprocw, global_my


    !%%% PACK: gg -> send_gg %%%
!$OMP parallel workshare
      send_gg(:,:,:,:,:) = (0._DP, 0._DP)
      send_psi(:,:,:,:) = (0._DP, 0._DP)
      send_chi(:,:,:,:) = (0._DP, 0._DP)
!$OMP end parallel workshare
!$OMP parallel do private(mx,my,iz,iv,im,izm,ibuff,iprocw)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                izm = (2*nz)*im + (iz + nz)
                ibuff = mod(izm, nbuff)
                iprocw = izm / nbuff
                send_gg(mx,my,iv,ibuff,iprocw) = gg(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do
      end do
!$OMP parallel do private(mx,my,iz,im,izm,ibuff,iprocw)
      do im = 0, nm
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                izm = (2*nz)*im + (iz + nz)
                ibuff = mod(izm, nbuff)
                iprocw = izm / nbuff
                send_psi(mx,my,ibuff,iprocw) = psi(mx,my,iz,im)
                send_chi(mx,my,ibuff,iprocw) = chi(mx,my,iz,im)
              end do
            end do
          end do
      end do

    !%%% TRANSPOSE: send_gg -> recv_gg %%%
      call MPI_Alltoall( send_gg,                      &
                         (2*nx+1)*(ny+1)*(2*nv)*nbuff, &
                         MPI_DOUBLE_COMPLEX,           &
                         recv_gg,                      &
                         (2*nx+1)*(ny+1)*(2*nv)*nbuff, &
                         MPI_DOUBLE_COMPLEX,           &
                         fft_comm_world,               &
                         ierr_mpi )
      call MPI_Alltoall( send_psi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         recv_psi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )
      call MPI_Alltoall( send_chi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         recv_chi,              &
                         (2*nx+1)*(ny+1)*nbuff, &
                         MPI_DOUBLE_COMPLEX,    &
                         fft_comm_world,        &
                         ierr_mpi )

    !%%% UNPACK: recv_gg -> wg %%%
!$OMP parallel do private(mx,global_my,iv,ibuff,iprocw,my)
      do ibuff = 0, nbuff-1
        do iv = 1, 2*nv
          do global_my = 0, global_ny
            do mx = -nx, nx
              iprocw = global_my / (ny+1)
              my = mod(global_my, ny+1)
              wg(mx,global_my,iv,ibuff) = recv_gg(mx,my,iv,ibuff,iprocw)
            end do
          end do
        end do
      end do
!$OMP parallel do private(mx,global_my,ibuff,iprocw,my)
      do ibuff = 0, nbuff-1
          do global_my = 0, global_ny
            do mx = -nx, nx
              iprocw = global_my / (ny+1)
              my = mod(global_my, ny+1)
              wps(mx,global_my,ibuff) = recv_psi(mx,my,ibuff,iprocw)
              wch(mx,global_my,ibuff) = recv_chi(mx,my,ibuff,iprocw)
            end do
          end do
      end do

                                         !%%% For debug %%%
                                         !iz = nz-1
                                         !im = nm
                                         !call MPI_Allgather(     &
                                         !     psi(-nx,0,iz,im),  &
                                         !     (2*nx+1)*(ny+1),   &
                                         !     MPI_DOUBLE_COMPLEX,&
                                         !     wch(-nx,0,0),      &
                                         !     (2*nx+1)*(ny+1),   &
                                         !     MPI_DOUBLE_COMPLEX,&
                                         !     fft_comm_world,    &
                                         !     ierr_mpi)
                                         !izm = (2*nz)*im+(iz+nz)
                                         !ibuff = mod(izm,nbuff)
                                         !iprocw = izm/nbuff
                                         !if (rankg==iprocw) then
                                         !write(888,*)"#",ibuff,iprocw
                                         !do my = 0, global_ny
                                         !  do mx = -nx, nx
                                         !    write(888,*) mx, my,    &
                                         !     dble(wps(mx,my,ibuff)),&
                                         !    aimag(wps(mx,my,ibuff)),&
                                         !     dble(wch(mx,my,0)),    &
                                         !    aimag(wch(mx,my,0))
                                         !  end do
                                         !  write(888,*)
                                         !end do
                                         !end if
                                         !call MPI_Finalize(ierr_mpi)
                                         !stop
                                         !%%%%%%%%%%%%%%%%%

  END SUBROUTINE trans_triad_transpose


!--------------------------------------
  SUBROUTINE trans_triad_coupling ( it, wg, wps, wch, jkpq_es, jpqk_es, jqkp_es,&
                                                      jkpq_em, jpqk_em, jqkp_em )
!--------------------------------------

    integer, intent(in)                              :: it
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny,1:2*nv,0:nbuff-1) :: wg
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:global_ny,0:nbuff-1)        :: wps, wch
    real(kind=DP), intent(out), &
      dimension(-nx:nx,-global_ny:global_ny) :: jkpq_es, jpqk_es, jqkp_es, &
                                                jkpq_em, jpqk_em, jqkp_em
    complex(kind=DP), dimension(-nx:nx,-global_ny:global_ny) :: gg, psi, chi
    real(kind=DP), dimension(-nx:nx,-global_ny:global_ny) :: &
                        wkpq_es, wpqk_es, wqkp_es, wkpq_em, wpqk_em, wqkp_em
    real(kind=DP) :: dky, wky(-global_ny:global_ny), cs1
    integer :: mx, my, px, py, qx, qy ! (mx,my) + (px,py) + (qx,qy) = (0,0)
    integer :: iv, ibuff

      !-set (mx,my)-
      mx = triad_diag_mxt(it)
      my = triad_diag_myt(it)

      !-set wky-
      dky = ky(1) - ky(0)
      do py = -global_ny, global_ny
        wky(py) = dky * real( py, kind=DP )
      end do

!$OMP parallel workshare
      wkpq_es(:,:) = 0._DP
      wpqk_es(:,:) = 0._DP
      wqkp_es(:,:) = 0._DP
      wkpq_em(:,:) = 0._DP
      wpqk_em(:,:) = 0._DP
      wqkp_em(:,:) = 0._DP
!$OMP end parallel workshare
      cs1 = sqrt( tau(ranks) / Anum(ranks) )
!$OMP parallel default(none) &
!$OMP shared(mx,my,kx,wky,vl,cs1,nbuff,wps,wch,wg) &
!$OMP shared(wkpq_es,wpqk_es,wqkp_es,wkpq_em,wpqk_em,wqkp_em) &
!$OMP private(px,py,qx,qy,iv,ibuff,psi,chi,gg)
!$OMP do reduction(+:wkpq_es) reduction(+:wpqk_es) reduction(+:wqkp_es) &
!$OMP    reduction(+:wkpq_em) reduction(+:wpqk_em) reduction(+:wqkp_em)
      do ibuff = 0, nbuff-1

        !-copy-
        do py = 0, global_ny
          do px = -nx, nx
            psi(px,py) = wps(px,py,ibuff)
            chi(px,py) = wch(px,py,ibuff)
          end do
        end do
        do py = 1, global_ny
          do px = -nx, nx
            psi(-px,-py) = conjg( wps(px,py,ibuff) )
            chi(-px,-py) = conjg( wch(px,py,ibuff) )
          end do
        end do

        do iv = 1, 2*nv

          !-copy-
          do py = 0, global_ny
            do px = -nx, nx
              gg(px,py) = wg(px,py,iv,ibuff)
            end do
          end do
          do py = 1, global_ny
            do px = -nx, nx
              gg(-px,-py) = conjg( wg(px,py,iv,ibuff) )
            end do
          end do

          !-triad coupling among (mx,my)+(px,py)+(qx,qy)=0-
          do py = max(-global_ny-my,-global_ny), min(global_ny,global_ny-my)
            qy = - my - py
            do px = max(-nx-mx,-nx), min(nx,nx-mx)
              qx = - mx - px
             !wkpq(px,py) = wkpq(px,py)                                         &
             !  - (- kx(px) * wky(qy) + wky(py) * kx(qx))                       &
             !  * real((  (psi(px,py) - cs1 * vl(iv) * chi(px,py)) * gg(qx,qy)  &
             !          - (psi(qx,qy) - cs1 * vl(iv) * chi(qx,qy)) * gg(px,py)) &
             !                                              * gg(mx,my), kind=DP)
             !wpqk(px,py) = wpqk(px,py)                                         &
             !  - (- kx(qx) * wky(my) + wky(qy) * kx(mx))                       &
             !  * real((  (psi(qx,qy) - cs1 * vl(iv) * chi(qx,qy)) * gg(mx,my)  &
             !          - (psi(mx,my) - cs1 * vl(iv) * chi(mx,my)) * gg(qx,qy)) &
             !                                              * gg(px,py), kind=DP)
             !wqkp(px,py) = wqkp(px,py)                                         &
             !  - (- kx(mx) * wky(py) + wky(my) * kx(px))                       &
             !  * real((  (psi(mx,my) - cs1 * vl(iv) * chi(mx,my)) * gg(px,py)  &
             !          - (psi(px,py) - cs1 * vl(iv) * chi(px,py)) * gg(mx,my)) &
             !                                              * gg(qx,qy), kind=DP)
              wkpq_es(px,py) = wkpq_es(px,py)             &
                - (- kx(px) * wky(qy) + wky(py) * kx(qx)) &
                * real((  (psi(px,py)) * gg(qx,qy)        &
                        - (psi(qx,qy)) * gg(px,py))       &
                                      * gg(mx,my), kind=DP)
              wpqk_es(px,py) = wpqk_es(px,py)             &
                - (- kx(qx) * wky(my) + wky(qy) * kx(mx)) &
                * real((  (psi(qx,qy)) * gg(mx,my)        &
                        - (psi(mx,my)) * gg(qx,qy))       &
                                      * gg(px,py), kind=DP)
              wqkp_es(px,py) = wqkp_es(px,py)             &
                - (- kx(mx) * wky(py) + wky(my) * kx(px)) &
                * real((  (psi(mx,my)) * gg(px,py)        &
                        - (psi(px,py)) * gg(mx,my))       &
                                      * gg(qx,qy), kind=DP)
              wkpq_em(px,py) = wkpq_em(px,py)                        &
                - (- kx(px) * wky(qy) + wky(py) * kx(qx))            &
                * real((  (- cs1 * vl(iv) * chi(px,py)) * gg(qx,qy)  &
                        - (- cs1 * vl(iv) * chi(qx,qy)) * gg(px,py)) &
                                                 * gg(mx,my), kind=DP)
              wpqk_em(px,py) = wpqk_em(px,py)                        &
                - (- kx(qx) * wky(my) + wky(qy) * kx(mx))            &
                * real((  (- cs1 * vl(iv) * chi(qx,qy)) * gg(mx,my)  &
                        - (- cs1 * vl(iv) * chi(mx,my)) * gg(qx,qy)) &
                                                 * gg(px,py), kind=DP)
              wqkp_em(px,py) = wqkp_em(px,py)                        &
                - (- kx(mx) * wky(py) + wky(my) * kx(px))            &
                * real((  (- cs1 * vl(iv) * chi(mx,my)) * gg(px,py)  &
                        - (- cs1 * vl(iv) * chi(px,py)) * gg(mx,my)) &
                                                 * gg(qx,qy), kind=DP)
            end do
          end do

        end do
      end do
!$OMP end do
!$OMP end parallel

      !-zz,vl,mu integration-
!$OMP parallel workshare
      jkpq_es(:,:) = 0._DP
      jpqk_es(:,:) = 0._DP
      jqkp_es(:,:) = 0._DP
      jkpq_em(:,:) = 0._DP
      jpqk_em(:,:) = 0._DP
      jqkp_em(:,:) = 0._DP
!$OMP end parallel workshare
      call MPI_Allreduce( wkpq_es, jkpq_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wpqk_es, jpqk_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wqkp_es, jqkp_es, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wkpq_em, jkpq_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wpqk_em, jpqk_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )
      call MPI_Allreduce( wqkp_em, jqkp_em, (2*nx+1)*(2*global_ny+1), &
                          MPI_DOUBLE_PRECISION,                       &
                          MPI_SUM, sub_comm_world, ierr_mpi )


  END SUBROUTINE trans_triad_coupling


END MODULE GKV_trans
