MODULE GKV_exb
!-------------------------------------------------------------------------------
!
!    E x B term
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fft,   only: fft_pre,  &
           fft_backward_Xfft, fft_backward_chXY, fft_backward_Yfft, &
           fft_forward_Yfft, fft_forward_chYX, fft_forward_Xfft
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

  public   exb_NL_term_es, exb_NL_term_em, exw_xw, eyw_xw, bxw_xw, byw_xw


CONTAINS


!--------------------------------------
  SUBROUTINE exb_NL_term_es( hh, psi, chi, ef )
!--------------------------------------
!  ExB nonlinear term calculation 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)                :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)                 :: ef

    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz)*(2*nv),0:nprocw-1,0:nm) :: send_df1, recv_df1
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz)*(2*nv),0:nprocw-1,0:nm) :: send_df2, recv_df2
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_exw, recv_exw
    complex(kind=DP),  &
      dimension(0:ny,0:nxw_size,(2*nz),0:nprocw-1,0:nm) :: send_eyw, recv_eyw
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

        do my = ist_y, iend_y
          do mx = 0, nx
            uikx(mx,my) = kx(mx) * ui
            uiky(mx,my) = ky(my) * ui
          end do
        end do

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
!$OMP shared (hh,psi,chi,ef)                          &
!$OMP shared (exbdf1,send_df1,recv_df1,exbdf1_xw) &
!$OMP shared (exbdf2,send_df2,recv_df2,exbdf2_xw) &
!$OMP shared (   exw,send_exw,recv_exw,   exw_xw) &
!$OMP shared (   eyw,send_eyw,recv_eyw,   eyw_xw) &
!$OMP shared (   bxw,                     bxw_xw) &
!$OMP shared (   byw,                     byw_xw) &
!$OMP shared (num_trans_exbdf,num_trans_exyw,num_trans_bxyw)     &
!$OMP private (im)

! -------------------------------------------------------------------
!     Backward FFT region
!     5 routines are overlaped:
!       exb_input, backward FFT in X, MPI_AlltoAll from X to Y,  
!       backward FFT in Y, real space cal.
! -------------------------------------------------------------------
      call exb_input (     hh(:,:,:,:,0), psi(:,:,:,0), chi(:,:,:,0),  &
                       exbdf1(:,:,:,:,0), exbdf2(:,:,:,:,0),           &
                            exw(:,:,:,0),      eyw(:,:,:,0),           &
                            bxw(:,:,:,0),      byw(:,:,:,0)           )
!$OMP barrier

      call exb_input (     hh(:,:,:,:,1), psi(:,:,:,1), chi(:,:,:,1),  &
                       exbdf1(:,:,:,:,1), exbdf2(:,:,:,:,1),           &
                            exw(:,:,:,1),      eyw(:,:,:,1),           &
                            bxw(:,:,:,1),      byw(:,:,:,1)           )

      call fft_backward_Xfft ( exbdf1(:,:,:,:,0), send_df1(:,:,:,:,0), num_trans_exbdf )
      call fft_backward_Xfft ( exbdf2(:,:,:,:,0), send_df2(:,:,:,:,0), num_trans_exbdf )
      call fft_backward_Xfft (      exw(:,:,:,0), send_exw(:,:,:,:,0), num_trans_exyw  )
      call fft_backward_Xfft (      eyw(:,:,:,0), send_eyw(:,:,:,:,0), num_trans_exyw  )
!$OMP barrier

      im = 0
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
!$OMP barrier

      im = 1
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )
!$OMP barrier

      do im = 2, nm-2
!$OMP master
        call fft_backward_chXY ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_df2(:,:,:,:,im), recv_df2(:,:,:,:,im), num_trans_exbdf )
        call fft_backward_chXY ( send_exw(:,:,:,:,im), recv_exw(:,:,:,:,im), num_trans_exyw  )
        call fft_backward_chXY ( send_eyw(:,:,:,:,im), recv_eyw(:,:,:,:,im), num_trans_exyw  )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
                         exbdf1(:,:,:,:,im+2), exbdf2(:,:,:,:,im+2),           &
                              exw(:,:,:,im+2),      eyw(:,:,:,im+2),           &
                              bxw(:,:,:,im+2),      byw(:,:,:,im+2)           )

        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
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
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_backward_Xfft ( exbdf1(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft ( exbdf2(:,:,:,:,im+1), send_df2(:,:,:,:,im+1), num_trans_exbdf )
        call fft_backward_Xfft (      exw(:,:,:,im+1), send_exw(:,:,:,:,im+1), num_trans_exyw  )
        call fft_backward_Xfft (      eyw(:,:,:,im+1), send_eyw(:,:,:,:,im+1), num_trans_exyw  )
  
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
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
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_backward_Yfft ( recv_df1(:,:,:,:,im-1), exbdf1_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,im-1), exbdf2_xw(:,:,:,:,im-1), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,im-1),      exw_xw(:,:,:,im-1), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,im-1),      eyw_xw(:,:,:,im-1), num_trans_exyw  )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
                                  exw_xw(:,:,:,im-2),      eyw_xw(:,:,:,im-2), & 
                                  bxw_xw(:,:,:,im-2),      byw_xw(:,:,:,im-2) )

        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,0), send_df1(:,:,:,:,0), num_trans_exbdf )
!$OMP barrier

      im = 0
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )

        call fft_backward_Yfft ( recv_df1(:,:,:,:,nm), exbdf1_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,nm),      exw_xw(:,:,:,nm), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,nm),      eyw_xw(:,:,:,nm), num_trans_exyw  )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,nm-1), exbdf2_xw(:,:,:,:,nm-1), &
                                  exw_xw(:,:,:,nm-1),      eyw_xw(:,:,:,nm-1), & 
                                  bxw_xw(:,:,:,nm-1),      byw_xw(:,:,:,nm-1) )
!$OMP barrier

      im = 1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), &
                                  exw_xw(:,:,:,nm),      eyw_xw(:,:,:,nm), &
                                  bxw_xw(:,:,:,nm),      byw_xw(:,:,:,nm) )
!$OMP barrier


! -------------------------------------------------------------------
!     Forward FFT region
!     4 routines are overlaped:
!       forward FFT in Y, MPI_AlltoAll form Y to X,
!       forward FFT in X, exb_output
! -------------------------------------------------------------------
      do im = 2, nm-1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
!$OMP barrier
      end do

      im = nm
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
!$OMP barrier

      call fft_forward_Xfft ( recv_df1(:,:,:,:,nm), exbdf1(:,:,:,:,nm), num_trans_exbdf )

      call exb_output ( exbdf1(:,:,:,:,nm-1), ef(:,:,:,:,nm-1) )
!$OMP barrier

      call exb_output ( exbdf1(:,:,:,:,nm), ef(:,:,:,:,nm) )
!$OMP end parallel


  END SUBROUTINE exb_NL_term_es


!--------------------------------------
  SUBROUTINE exb_NL_term_em( hh, psi, chi, ef )
!--------------------------------------
!  ExB nonlinear term calculation 

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)                :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)                 :: ef

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

        do my = ist_y, iend_y
          do mx = 0, nx
            uikx(mx,my) = kx(mx) * ui
            uiky(mx,my) = ky(my) * ui
          end do
        end do

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
!$OMP shared (hh,psi,chi,ef)                          &
!$OMP shared (exbdf1,send_df1,recv_df1,exbdf1_xw) &
!$OMP shared (exbdf2,send_df2,recv_df2,exbdf2_xw) &
!$OMP shared (   exw,send_exw,recv_exw,   exw_xw) &
!$OMP shared (   eyw,send_eyw,recv_eyw,   eyw_xw) &
!$OMP shared (   bxw,send_bxw,recv_bxw,   bxw_xw) &
!$OMP shared (   byw,send_byw,recv_byw,   byw_xw) &
!$OMP shared (num_trans_exbdf,num_trans_exyw,num_trans_bxyw)     &
!$OMP private (im)

! -------------------------------------------------------------------
!     Backward FFT region
!     5 routines are overlaped:
!       exb_input, backward FFT in X, MPI_AlltoAll from X to Y,  
!       backward FFT in Y, real space cal.
! -------------------------------------------------------------------
      call exb_input (     hh(:,:,:,:,0), psi(:,:,:,0), chi(:,:,:,0),  &
                       exbdf1(:,:,:,:,0), exbdf2(:,:,:,:,0),           &
                            exw(:,:,:,0),      eyw(:,:,:,0),           &
                            bxw(:,:,:,0),      byw(:,:,:,0)           )
!$OMP barrier

      call exb_input (     hh(:,:,:,:,1), psi(:,:,:,1), chi(:,:,:,1),  &
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
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
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
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
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
        call exb_input (     hh(:,:,:,:,im+2), psi(:,:,:,im+2), chi(:,:,:,im+2),  &
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

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
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

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
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

        call exb_realspcal ( exbdf1_xw(:,:,:,:,im-2), exbdf2_xw(:,:,:,:,im-2), &
                                  exw_xw(:,:,:,im-2),      eyw_xw(:,:,:,im-2), & 
                                  bxw_xw(:,:,:,im-2),      byw_xw(:,:,:,im-2) )

        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,0), send_df1(:,:,:,:,0), num_trans_exbdf )
!$OMP barrier

      im = 0
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )

        call fft_backward_Yfft ( recv_df1(:,:,:,:,nm), exbdf1_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_df2(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), num_trans_exbdf )
        call fft_backward_Yfft ( recv_exw(:,:,:,:,nm),      exw_xw(:,:,:,nm), num_trans_exyw  )
        call fft_backward_Yfft ( recv_eyw(:,:,:,:,nm),      eyw_xw(:,:,:,nm), num_trans_exyw  )
        call fft_backward_Yfft ( recv_bxw(:,:,:,:,nm),      bxw_xw(:,:,:,nm), num_trans_bxyw  )
        call fft_backward_Yfft ( recv_byw(:,:,:,:,nm),      byw_xw(:,:,:,nm), num_trans_bxyw  )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,nm-1), exbdf2_xw(:,:,:,:,nm-1), &
                                  exw_xw(:,:,:,nm-1),      eyw_xw(:,:,:,nm-1), & 
                                  bxw_xw(:,:,:,nm-1),      byw_xw(:,:,:,nm-1) )
!$OMP barrier

      im = 1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_realspcal ( exbdf1_xw(:,:,:,:,nm), exbdf2_xw(:,:,:,:,nm), &
                                  exw_xw(:,:,:,nm),      eyw_xw(:,:,:,nm), &
                                  bxw_xw(:,:,:,nm),      byw_xw(:,:,:,nm) )
!$OMP barrier


! -------------------------------------------------------------------
!     Forward FFT region
!     4 routines are overlaped:
!       forward FFT in Y, MPI_AlltoAll form Y to X,
!       forward FFT in X, exb_output
! -------------------------------------------------------------------
      do im = 2, nm-1
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Yfft ( exbdf1_xw(:,:,:,:,im+1), send_df1(:,:,:,:,im+1), num_trans_exbdf )
  
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
!$OMP barrier
      end do

      im = nm
!$OMP master
        call fft_forward_chYX ( send_df1(:,:,:,:,im), recv_df1(:,:,:,:,im), num_trans_exbdf )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call fft_forward_Xfft ( recv_df1(:,:,:,:,im-1), exbdf1(:,:,:,:,im-1), num_trans_exbdf )

        call exb_output ( exbdf1(:,:,:,:,im-2), ef(:,:,:,:,im-2) )
!$OMP barrier

      call fft_forward_Xfft ( recv_df1(:,:,:,:,nm), exbdf1(:,:,:,:,nm), num_trans_exbdf )

      call exb_output ( exbdf1(:,:,:,:,nm-1), ef(:,:,:,:,nm-1) )
!$OMP barrier

      call exb_output ( exbdf1(:,:,:,:,nm), ef(:,:,:,:,nm) )
!$OMP end parallel


  END SUBROUTINE exb_NL_term_em


!--------------------------------------
  SUBROUTINE exb_input ( hh, psi, chi, wkdf1, wkdf2, wkexw, wkeyw, wkbxw, wkbyw )
!--------------------------------------
!     Data input for E x B term calculation

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)    :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb)   :: psi, chi
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

      do iv = 1, 2*nv
!$OMP do schedule (dynamic)
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
!$OMP end do nowait
      end do


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


  END SUBROUTINE exb_input


!--------------------------------------
  SUBROUTINE exb_realspcal ( wkdf1_xw, wkdf2_xw, wkexw_xw, wkeyw_xw, wkbxw_xw, wkbyw_xw )
!--------------------------------------
!     Calculate E x B term in real space

    complex(kind=DP), intent(inout), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1,1:2*nv) :: wkdf1_xw
    complex(kind=DP), intent(in), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1,1:2*nv) :: wkdf2_xw
    complex(kind=DP), intent(in), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1)        :: wkexw_xw, wkeyw_xw
    complex(kind=DP), intent(in), &
      dimension(0:nyw,0:nxw_size,-nz:nz-1)        :: wkbxw_xw, wkbyw_xw
    real(kind=DP) ::  cef, cs1
    integer  ::  mx, my, iz, iv

      cef = 1._DP / real( nnx*nny, kind=DP )
      cs1 = sqrt( tau(ranks) / Anum(ranks) )

! --- Real space calculation
!$OMP master
                                           call clock_sta(1430)
                                         ! call fapp_start("nlterm_realspcal",1430,1)
!$OMP end master
      do iv = 1, 2*nv
!$OMP do schedule (dynamic)
        do iz = -nz, nz-1
          do mx = ist_xw, iend_xw
            do my = 0, nyw
              wkdf1_xw(my,mx,iz,iv) =                           &
                 cmplx(      real ( wkdf1_xw(my,mx,iz,iv), kind=DP ) &
                           * real ( wkeyw_xw(my,mx,iz) - cs1 * vl(iv) &
                                  * wkbyw_xw(my,mx,iz),    kind=DP ) &
                       -     real ( wkdf2_xw(my,mx,iz,iv), kind=DP ) &
                           * real ( wkexw_xw(my,mx,iz) - cs1 * vl(iv) &
                                  * wkbxw_xw(my,mx,iz),    kind=DP ) &
                       ,    aimag ( wkdf1_xw(my,mx,iz,iv)          ) &
                          * aimag ( wkeyw_xw(my,mx,iz) - cs1 * vl(iv) &
                                  * wkbyw_xw(my,mx,iz)             ) &
                      -     aimag ( wkdf2_xw(my,mx,iz,iv)          ) &
                          * aimag ( wkexw_xw(my,mx,iz) - cs1 * vl(iv) &
                                  * wkbxw_xw(my,mx,iz)             ) &
                      , kind=DP ) * cef
            end do
          end do
        end do
!$OMP end do nowait
      end do
!$OMP master
                                         ! call fapp_stop("nlterm_realspcal",1430,1)
                                           call clock_end(1430)
!$OMP end master


  END SUBROUTINE exb_realspcal


!--------------------------------------
  SUBROUTINE exb_output ( wkdf1, wf )
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
      do iv = 1, 2*nv
!$OMP do schedule (dynamic)
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
!$OMP end do nowait
      end do
!$OMP master
                                         ! call fapp_stop("nlterm_output",1450,1)
                                           call clock_end(1450)
!$OMP end master


  END SUBROUTINE exb_output


END MODULE GKV_exb
