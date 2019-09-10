MODULE GKV_advnc
!-------------------------------------------------------------------------------
!
!    Flux surface and field line averages
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_fld,   only: fld_esfield, fld_emfield_hh, fld_hh2ff
  use GKV_exb,   only: exb_NL_term_es, exb_NL_term_em
  use GKV_colli, only: colli_LB_model, colli_GK_CT, colli_GK_CT6, colli_GK_CF_DT, & 
                       colli_moment_calc, colli_moment_redc, colli_comm_alltoall, & 
                       colli_dfdvp6, colli_dfdvp, colli_zeroset, &
                       colli_hhset, colli_wwset
  use GKV_bndry, only: bndry_bound_e,  &
      bndry_bound_f_buffin, bndry_bound_f_sendrecv, bndry_bound_f_buffout,  &
      bndry_shifts_v_buffin, bndry_shifts_v_sendrecv, bndry_shifts_v_buffout,  &
      bndry_zv_sendrecv, bndry_zvm_bound_f, &
      bndry_shifts_m_buffin, bndry_shifts_m_sendrecv, bndry_shifts_m_buffout
  use GKV_clock, only: clock_sta, clock_end
  use GKV_zfilter, only: zfilter
  use GKV_tips,  only: tips_reality

  implicit none

  private

  public   advnc_rkgsteps, caldlt_LB, caldlt_full


CONTAINS


!--------------------------------------
  SUBROUTINE advnc_rkgsteps( ff, phi, Al, hh, dh, cft, cff, ef )
!--------------------------------------
!     time integration of GK equation using RKG method

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh, dh, cft, cff, ef

    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: qh
    integer :: istep, mx, my, iz, iv, im, ierr_mpi


!$OMP parallel do collapse(2)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  qh(mx,my,iz,iv,im) = ( 0._DP, 0._DP )
                end do
              end do
            end do
          end do
        end do

      do istep = 1, 4

                                           call clock_sta(11)
                                         ! call fapp_start("rkg",11,1)
        call rkg( hh, dh, qh, istep )
                                         ! call fapp_stop("rkg",11,1)
                                           call clock_end(11)

        call tips_reality ( hh )

                                           call clock_sta(12)
                                         ! call fapp_start("esfield",12,1)
        if ( beta .ne. 0._DP ) then
          call fld_emfield_hh( hh, Al )
        end if
        call fld_hh2ff( hh, Al, ff )
        call fld_esfield( ff, phi )
                                         ! call fapp_stop("esfield",12,1)
                                           call clock_end(12)

      if ( trim(col_type) == "LB" ) then
        call caldlt_LB( ff, phi, Al, hh, dh, cft, cff, ef )
      else if ( trim(col_type)  == "full" ) then
        call caldlt_full( ff, phi, Al, hh, dh, cft, cff, ef )
      else 
        write(olog,*) "## Illegal choice for col_type!! ---> stop"
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop
      end if

      end do


  END SUBROUTINE advnc_rkgsteps


!--------------------------------------
  SUBROUTINE rkg( hh, dh, qh, istep )
!--------------------------------------
!     Runge-Kutta-Gill

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: dh, qh

    integer, intent(in) :: istep
    real(kind=DP) :: c1, c2, cq, c0
    integer :: mx, my, iz, iv, im


      if      ( istep == 1 ) then
        c1   =  0.5_DP
        c2   = -1._DP
        cq   = -2._DP
        c0   =  1._DP
      else if ( istep == 2 ) then
        c1   =  1._DP - sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 3 ) then
        c1   =  1._DP + sqrt( 0.5_DP )
        c2   = -c1
        cq   =  1._DP - 3._DP * c1
        c0   =  2._DP * c1
      else if ( istep == 4 ) then
        c1   =  1._DP / 6._DP
        c2   = -1._DP / 3._DP
        cq   =  0._DP
        c0   =  0._DP
      end if


!$OMP parallel do collapse(2)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                hh(mx,my,iz,iv,im)   =      hh(mx,my,iz,iv,im) &
                                     + c1 * dt * dh(mx,my,iz,iv,im) &
                                     + c2 * qh(mx,my,iz,iv,im)
                qh(mx,my,iz,iv,im)   = cq * qh(mx,my,iz,iv,im) &
                                     + c0 * dt * dh(mx,my,iz,iv,im)
              end do
            end do
          end do
        end do
      end do


  END SUBROUTINE rkg


!--------------------------------------
  SUBROUTINE caldlt_LB( ff, phi, Al, hh, dh, cft, cff, ef )
!--------------------------------------
!     increment of delta-f within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: dh, cff, cft, ef

    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: lf
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb1_bottom, zb1_top
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb2_bottom, zb2_top
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm)   :: vb1, vb2
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1, mb2
    integer  ::  mx, my, iz, iv, im


                                           call clock_sta(13)
                                         ! call fapp_start("literm",13,1)
                                           call clock_sta(1340)
                                         ! call fapp_start("literm_other",1340,1)
!$OMP parallel do collapse(2)
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
                                         ! call fapp_stop("literm_other",1340,1)
                                           call clock_end(1340)

      call bndry_bound_e ( psi )
                                         ! call fapp_stop("literm",13,1)
                                           call clock_end(13)

                                           call clock_sta(14)
                                         ! call fapp_start("nlterm",14,1)
      if( trim(calc_type) == "nonlinear" ) then
        if ( beta == 0._DP ) then
          call exb_NL_term_es( hh, psi, chi, ef )
        else
          call exb_NL_term_em( hh, psi, chi, ef )
        end if
      else
!$OMP parallel do collapse(2)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ef(mx,my,iz,iv,im) = ( 0._DP, 0._DP )
                end do
              end do
            end do
          end do
        end do
      end if
                                         ! call fapp_stop("nlterm",14,1)
                                           call clock_end(14)

! --- Annihilation test for maxwellian
      if (icheck == 1) then
!$OMP parallel do collapse(2)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
!                  ff(mx,my,iz,iv,im) = fmx(iz,iv,im)
                  ff(mx,my,iz,iv,im) = dexp( -xxa(iz,iv,im)**2 )
                end do
              end do
            end do
          end do
        end do
!        call bndry_zvm_bound_f( ff )
      end if 
! --- 
                                           call clock_sta(13)
                                         ! call fapp_start("literm",13,1)
!$OMP parallel default (none) &
!$OMP shared(ff,psi,chi,lf,dh,ef) &
!$OMP shared(cft,cff) &
!$OMP shared(zb1_bottom,zb1_top,zb2_bottom,zb2_top,vb1,vb2,mb1,mb2) &
!$OMP shared(rankw,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im)
      call bndry_shifts_m_buffin ( ff, mb1, mb2 )
!$OMP barrier

!$OMP master
      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!$OMP end master
!!!!! barrier ! for test without overlaps
      call literm_k ( ff, psi, chi, lf )
      do im = 0, nm
        call bndry_bound_f_buffin ( ff(:,:,:,:,im), zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im) )
        call bndry_shifts_v_buffin ( ff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
      end do
      call colli_zeroset( cff )
!$OMP barrier

      im = 0
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call bndry_shifts_m_buffout ( mb2, ff )
!$OMP barrier

      im = 1
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1), &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!$OMP barrier

      im = 2
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1),  &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
        call colli_LB_model( ff, im-2, cft(:,:,:,:,im-2) )
        call literm_zv ( ff(:,:,:,:,im-2), psi(:,:,:,im-2), im-2, lf(:,:,:,:,im-2) )
!$OMP barrier

      do im = 3, nm
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps
        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1),  &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
        call colli_LB_model( ff, im-2, cft(:,:,:,:,im-2) )
        call literm_zv ( ff(:,:,:,:,im-2), psi(:,:,:,im-2), im-2, lf(:,:,:,:,im-2) )
        call calc_dh ( lf(:,:,:,:,im-3), ef(:,:,:,:,im-3),  &
                       cft(:,:,:,:,im-3), cff(:,:,:,:,im-3), dh(:,:,:,:,im-3) )
!$OMP barrier
      end do

      call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,nm), zb2_top(:,:,:,:,nm),  &
                                   ff(:,:,:,:,nm) )
      call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), ff(:,:,:,:,nm) )
      call colli_LB_model( ff, nm-1, cft(:,:,:,:,nm-1) )
      call literm_zv ( ff(:,:,:,:,nm-1), psi(:,:,:,nm-1), nm-1, lf(:,:,:,:,nm-1) )
      call calc_dh ( lf(:,:,:,:,nm-2), ef(:,:,:,:,nm-2),  &
                     cft(:,:,:,:,nm-2), cff(:,:,:,:,nm-2), dh(:,:,:,:,nm-2) )
!$OMP barrier

      call colli_LB_model( ff, nm, cft(:,:,:,:,nm) )
      call literm_zv ( ff(:,:,:,:,nm), psi(:,:,:,nm), nm, lf(:,:,:,:,nm) )
      call calc_dh ( lf(:,:,:,:,nm-1), ef(:,:,:,:,nm-1),  &
                     cft(:,:,:,:,nm-1), cff(:,:,:,:,nm-1), dh(:,:,:,:,nm-1) )
!$OMP barrier

      call calc_dh ( lf(:,:,:,:,nm), ef(:,:,:,:,nm),  &
                     cft(:,:,:,:,nm), cff(:,:,:,:,nm), dh(:,:,:,:,nm) )
!$OMP end parallel
                                         ! call fapp_stop("literm",13,1)
                                           call clock_end(13)


                                           call clock_sta(15)
                                         ! call fapp_start("zfilter",15,1)
      if ( trim(z_filt) == "on" ) then
        call zfilter ( dh )
      end if                                
                                         ! call fapp_stop("zfilter",15,1)
                                           call clock_end(15)


  END SUBROUTINE caldlt_LB


!--------------------------------------
  SUBROUTINE caldlt_full( ff, phi, Al, hh, dh, cft, cff, ef )
!--------------------------------------
!     increment of delta-f within a time step

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi, Al
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: hh
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: dh, cff, cft, ef

    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: lf
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm) :: psi, chi
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb1_bottom, zb1_top
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,0:nzb-1,1:2*nv,0:nm) :: zb2_bottom, zb2_top
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nvb,0:nm)   :: vb1, vb2
    complex(kind=DP),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,1:2*nvb) :: mb1, mb2

    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6) :: wrkm, moment_ab, moment_ba
    complex(kind=DP), &
      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1) :: moment_ab_wk, moment_ba_wk
    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: dfdvp

    integer  ::  nn
    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh2
    integer  ::  mx, my, iz, iv, im, ii


                                           call clock_sta(13)
                                         ! call fapp_start("literm",13,1)
                                           call clock_sta(1340)
                                         ! call fapp_start("literm_other",1340,1)
!$OMP parallel do collapse(2)
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
                                         ! call fapp_stop("literm_other",1340,1)
                                           call clock_end(1340)

      call bndry_bound_e ( psi )
                                         ! call fapp_stop("literm",13,1)
                                           call clock_end(13)

                                           call clock_sta(14)
                                         ! call fapp_start("nlterm",14,1)
      if( trim(calc_type) == "nonlinear" ) then
        if ( beta == 0._DP ) then
          call exb_NL_term_es( hh, psi, chi, ef )
        else
          call exb_NL_term_em( hh, psi, chi, ef )
        end if
      else
!$OMP parallel do collapse(2)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ef(mx,my,iz,iv,im) = ( 0._DP, 0._DP )
                end do
              end do
            end do
          end do
        end do
      end if
                                         ! call fapp_stop("nlterm",14,1)
                                           call clock_end(14)

! --- Annihilation test for maxwellian
      if (icheck == 1) then
!$OMP parallel do collapse(2)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ff(mx,my,iz,iv,im) = fmx(iz,iv,im)
!                  ff(mx,my,iz,iv,im) = dexp( -xxa(iz,iv,im)**2 )
                end do
              end do
            end do
          end do
        end do
      end if 
! --- 

                                           call clock_sta(13)
                                         ! call fapp_start("literm",13,1)
!$OMP parallel default (none) &
!$OMP shared(ff,phi,dfdvp) &
!$OMP shared(mb1,mb2,vb1,vb2) &
!$OMP shared(wrkm,moment_ab,moment_ba,hh2) &
!$OMP shared(cft,cff) &
!$OMP private(im,ii) 

      call bndry_shifts_m_buffin ( ff, mb1, mb2 )
      call colli_hhset(hh2,phi,ff)
      call colli_wwset(wrkm)
!$OMP barrier

!----------------------------------------------------- ovlp1
!$OMP master
      call bndry_shifts_m_sendrecv ( mb1, mb2 )
!$OMP end master
!!!!! barrier ! for test without overlaps
      call colli_moment_calc( hh2, phi, wrkm )
      call colli_zeroset( cff )

!$OMP barrier
!-----------------------------------------------------

!----------------------------------------------------- non-ovlp part
      call bndry_shifts_m_buffout ( mb2, ff )
!$OMP barrier
!-----------------------------------------------------

!----------------------------------------------------- ovlp2
!$OMP master
      call colli_moment_redc( wrkm, moment_ab )
!$OMP end master
!!!!! barrier ! for test without overlaps

      call colli_dfdvp6( ff, dfdvp )  ! 6th-order CFD
!!!      call colli_dfdvp( ff, dfdvp )  ! 4th-order CFD
!$OMP barrier
!-----------------------------------------------------

!----------------------------------------------------- non-ovlp part
      do im = 0, nm
        call bndry_shifts_v_buffin ( dfdvp(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
      end do
!$OMP barrier
!-----------------------------------------------------

!----------------------------------------------------- ovlp3
      im = 0
!$OMP master
        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

 ! --- No calculations appear here in f0.48 (Nakata July2015)
 !!!      call colli_zeroset( cdt )

!$OMP barrier

      do im = 1, nm
!$OMP master
        call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), dfdvp(:,:,:,:,im-1) )
!$OMP barrier
      end do

        call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), dfdvp(:,:,:,:,nm) )

!-----------------------------------------------------
!$OMP end parallel


!$OMP parallel default (none) &
!$OMP shared(ff,psi,chi,lf,dh,ef) &
!$OMP shared(cft,cff) &
!$OMP shared(phi,dfdvp,moment_ab,moment_ba,wrkm,moment_ab_wk,moment_ba_wk) &
!$OMP shared(zb1_bottom,zb1_top,zb2_bottom,zb2_top,vb1,vb2,mb1,mb2) &
!$OMP shared(rankw,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,nn) 

!----------------------------------------------------- ovlp4
!$OMP master
      call colli_comm_alltoall( moment_ab, moment_ba )
!$OMP end master
!!!!! barrier ! for test without overlaps

      do im = 0, nm
        call bndry_bound_f_buffin ( ff(:,:,:,:,im), zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im) )
        call bndry_shifts_v_buffin ( ff(:,:,:,:,im), vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
      end do
!$OMP barrier

!$OMP do collapse(2)
do nn = 0, ns-1
  do iz = -nz, nz-1
    do my = ist_y, iend_y
      do mx = -nx, nx
        moment_ab_wk(1,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,1)
        moment_ab_wk(2,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,2)
        moment_ab_wk(3,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,3)
        moment_ab_wk(4,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,4)
        moment_ab_wk(5,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,5)
        moment_ab_wk(6,mx,my,iz,nn) = moment_ab(mx,my,iz,nn,6)

        moment_ba_wk(1,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,1)
        moment_ba_wk(2,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,2)
        moment_ba_wk(3,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,3)
        moment_ba_wk(4,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,4)
        moment_ba_wk(5,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,5)
        moment_ba_wk(6,mx,my,iz,nn) = moment_ba(mx,my,iz,nn,6)
      enddo
    enddo
  enddo
enddo
!$OMP enddo
!-----------------------------------------------------

!----------------------------------------------------- ovlp5
      im = 0
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

        call literm_k ( ff, psi, chi, lf )
!$OMP barrier
!----------------------------------------------------- 

!----------------------------------------------------- ovlp6
      im = 1
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1), &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
!$OMP barrier
!----------------------------------------------------- 

!----------------------------------------------------- ovlp7
      im = 2
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1),  &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!        call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD 
        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
        call literm_zv ( ff(:,:,:,:,im-2), psi(:,:,:,im-2), im-2, lf(:,:,:,:,im-2) )
!$OMP barrier
!----------------------------------------------------- 

!----------------------------------------------------- ovlp8
      do im = 3, nm
!$OMP master
      !  call bndry_bound_f_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
      !                                zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im) )
      !  call bndry_shifts_v_sendrecv ( vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
        call bndry_zv_sendrecv ( zb1_bottom(:,:,:,:,im), zb1_top(:,:,:,:,im),  &
                                 zb2_bottom(:,:,:,:,im), zb2_top(:,:,:,:,im),  &
                                 vb1(:,:,:,:,im), vb2(:,:,:,:,im) )
!$OMP end master
!!!!! barrier ! for test without overlaps

        call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,im-1), zb2_top(:,:,:,:,im-1),  &
                                     ff(:,:,:,:,im-1) )
        call bndry_shifts_v_buffout ( vb2(:,:,:,:,im-1), ff(:,:,:,:,im-1) )
        call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 6th-order CFD
!!!        call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,im-2), im-2, cft(:,:,:,:,im-2) )  ! 4th-order CFD
        call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, im-2, cff(:,:,:,:,im-2) )
        call literm_zv ( ff(:,:,:,:,im-2), psi(:,:,:,im-2), im-2, lf(:,:,:,:,im-2) )
        call calc_dh ( lf(:,:,:,:,im-3), ef(:,:,:,:,im-3),  &
                       cft(:,:,:,:,im-3), cff(:,:,:,:,im-3), dh(:,:,:,:,im-3) )
!$OMP barrier
      end do
!----------------------------------------------------- 

      call bndry_bound_f_buffout ( zb2_bottom(:,:,:,:,nm), zb2_top(:,:,:,:,nm),  &
                                   ff(:,:,:,:,nm) )
      call bndry_shifts_v_buffout ( vb2(:,:,:,:,nm), ff(:,:,:,:,nm) )
      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 6th-order CFD
!!!      call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,nm-1), nm-1, cft(:,:,:,:,nm-1) )  ! 4th-order CFD
      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm-1, cff(:,:,:,:,nm-1) )
      call literm_zv ( ff(:,:,:,:,nm-1), psi(:,:,:,nm-1), nm-1, lf(:,:,:,:,nm-1) )
      call calc_dh ( lf(:,:,:,:,nm-2), ef(:,:,:,:,nm-2),  &
                     cft(:,:,:,:,nm-2), cff(:,:,:,:,nm-2), dh(:,:,:,:,nm-2) )
!$OMP barrier

      call colli_GK_CT6( ff, phi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 6th-order CFD
!!!      call colli_GK_CT( ff, phi, dfdvp(:,:,:,:,nm), nm, cft(:,:,:,:,nm) )  ! 4th-order CFD
      call colli_GK_CF_DT( moment_ba_wk, moment_ab_wk, nm, cff(:,:,:,:,nm) )
      call literm_zv ( ff(:,:,:,:,nm), psi(:,:,:,nm), nm, lf(:,:,:,:,nm) )
      call calc_dh ( lf(:,:,:,:,nm-1), ef(:,:,:,:,nm-1),  &
                     cft(:,:,:,:,nm-1), cff(:,:,:,:,nm-1), dh(:,:,:,:,nm-1) )
!$OMP barrier

      call calc_dh ( lf(:,:,:,:,nm), ef(:,:,:,:,nm),  &
                     cft(:,:,:,:,nm), cff(:,:,:,:,nm), dh(:,:,:,:,nm) )

!$OMP end parallel
                                         ! call fapp_stop("literm",13,1)
                                           call clock_end(13)

                                           call clock_sta(15)
                                         ! call fapp_start("zfilter",15,1)
      if ( trim(z_filt) == "on" ) then
        call zfilter ( dh )
      end if                                
                                         ! call fapp_stop("zfilter",15,1)
                                           call clock_end(15)


  END SUBROUTINE caldlt_full


!--------------------------------------
  SUBROUTINE literm_k ( ff, psi, chi, lf )
!--------------------------------------
!     z-derivative of ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,0:nm)                :: psi, chi
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)           :: lf

    real(kind=DP) :: cs1, cs2
    integer  ::  mx, my, iz, iv, im


!$OMP master
                                           call clock_sta(1320)
                                         ! call fapp_start("literm_perp",1320,1)
!$OMP end master

      cs1    = sgn(ranks) * Znum(ranks) / tau(ranks)
      cs2    = sqrt( tau(ranks) / Anum(ranks) )

!$OMP do schedule (dynamic)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                lf(mx,my,iz,iv,im) =                                  &
                   - ui * kvd(mx,my,iz,iv,im) * ff(mx,my,iz,iv,im)    &
                   - cs1 * fmx(iz,iv,im) * (                          &
                       + ui * kvd(mx,my,iz,iv,im) * psi(mx,my,iz,im)  &
                       - ui * kvs(mx,my,iz,iv,im)                     &
                            * ( psi(mx,my,iz,im) - cs2 * vl(iv) * chi(mx,my,iz,im) ) )
              end do
            end do
          end do
        end do
      end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_perp",1320,1)
                                           call clock_end(1320)
!$OMP end master


  END SUBROUTINE literm_k


! --- CFD for ff and psi
!--------------------------------------
  SUBROUTINE literm_zv ( ff, psi, im, lf )
!--------------------------------------
!     z-derivative of ff

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb)            :: psi
    integer, intent(in) :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)             :: lf

    real(kind=DP), dimension(-nz:nz-1) :: cefz, cefz2
    real(kind=DP) :: cefv, cs1
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1330)
                                         ! call fapp_start("literm_para",1330,1)
!$OMP end master

      cs1    = sgn(ranks) * Znum(ranks) / tau(ranks)
      do iz = -nz, nz-1
        cefz(iz)   = 1._DP / ( 12._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
        cefz2(iz)  = 1._DP / ( 60._DP * dpara(iz) ) * sqrt( tau(ranks) / Anum(ranks) )
      end do
      cefv   = 1._DP / ( 12._DP * dv ) * sqrt( tau(ranks) / Anum(ranks) )

      if (trim(z_calc) == "cf4") then

!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                 - vl(iv) * cefz(iz) * (              &
                     -         ff(mx,my,iz+2,iv)      &
                     + 8._DP * ff(mx,my,iz+1,iv)      &
                     - 8._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )    &
                 + mir(iz,im) * cefv * (              &
                     -         ff(mx,my,iz,iv+2)      &
                     + 8._DP * ff(mx,my,iz,iv+1)      &
                     - 8._DP * ff(mx,my,iz,iv-1)      &
                     +         ff(mx,my,iz,iv-2) )    &
                 - cs1 * fmx(iz,iv,im) * (            &
                       vl(iv) * cefz(iz) * (          &
                         -         psi(mx,my,iz+2)    &
                         + 8._DP * psi(mx,my,iz+1)    &
                         - 8._DP * psi(mx,my,iz-1)    &
                         +         psi(mx,my,iz-2) ) )&
                 - art_diff * (                       &
                     +         ff(mx,my,iz+2,iv)      &
                     - 4._DP * ff(mx,my,iz+1,iv)      &
                     + 6._DP * ff(mx,my,iz  ,iv)      &
                     - 4._DP * ff(mx,my,iz-1,iv)      &
                     +         ff(mx,my,iz-2,iv) )
            end do
          end do
        end do
      end do
!$OMP end do nowait

      else if (trim(z_calc) == "up5") then

        do iv = 1, 2*nv
          if ( vl(iv) > 0._DP ) then
!$OMP do schedule (dynamic)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         - 3._DP * ff(mx,my,iz+2,iv)      &
                         +30._DP * ff(mx,my,iz+1,iv)      &
                         +20._DP * ff(mx,my,iz  ,iv)      &
                         -60._DP * ff(mx,my,iz-1,iv)      &
                         +15._DP * ff(mx,my,iz-2,iv)      &
                         - 2._DP * ff(mx,my,iz-3,iv) )    &
                     + mir(iz,im) * cefv * (              &
                         -         ff(mx,my,iz,iv+2)      &
                         + 8._DP * ff(mx,my,iz,iv+1)      &
                         - 8._DP * ff(mx,my,iz,iv-1)      &
                         +         ff(mx,my,iz,iv-2) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz(iz) * (          &
                             -         psi(mx,my,iz+2)    &
                             + 8._DP * psi(mx,my,iz+1)    &
                             - 8._DP * psi(mx,my,iz-1)    &
                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          else
!$OMP do schedule (dynamic)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  lf(mx,my,iz,iv) = lf(mx,my,iz,iv)       &
                     - vl(iv) * cefz2(iz) * (             &
                         + 2._DP * ff(mx,my,iz+3,iv)      &
                         -15._DP * ff(mx,my,iz+2,iv)      &
                         +60._DP * ff(mx,my,iz+1,iv)      &
                         -20._DP * ff(mx,my,iz  ,iv)      &
                         -30._DP * ff(mx,my,iz-1,iv)      &
                         + 3._DP * ff(mx,my,iz-2,iv) )    &
                     + mir(iz,im) * cefv * (              &
                         -         ff(mx,my,iz,iv+2)      &
                         + 8._DP * ff(mx,my,iz,iv+1)      &
                         - 8._DP * ff(mx,my,iz,iv-1)      &
                         +         ff(mx,my,iz,iv-2) )    &
                     - cs1 * fmx(iz,iv,im) * (            &
                           vl(iv) * cefz(iz) * (          &
                             -         psi(mx,my,iz+2)    &
                             + 8._DP * psi(mx,my,iz+1)    &
                             - 8._DP * psi(mx,my,iz-1)    &
                             +         psi(mx,my,iz-2) ) )
                end do
              end do
            end do
!$OMP end do nowait
          end if
        end do

      else

        write(olog,*) "## Illegal choice for z_calc!! ---> stop"
        call flush(olog)
        call MPI_Finalize(ierr_mpi)
        stop

      end if

!$OMP master
                                         ! call fapp_stop("literm_para",1330,1)
                                           call clock_end(1330)
!$OMP end master


  END SUBROUTINE literm_zv


!--------------------------------------
  SUBROUTINE calc_dh ( lf, ef, cft, cff, dh )
!--------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)             :: lf, ef, cft, cff
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)             :: dh

    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1340)
                                         ! call fapp_start("literm_other",1340,1)
!$OMP end master

!$OMP do schedule (dynamic)
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              dh(mx,my,iz,iv) = ( lf(mx,my,iz,iv)  &
                                - ef(mx,my,iz,iv)  &
                                + cft(mx,my,iz,iv) &
                                + cff(mx,my,iz,iv) )
            end do
          end do
          if( rankw == 0 )  then
              dh(0,0,iz,iv)   = ( 0._DP, 0._DP )
          end if
        end do
      end do
!$OMP end do nowait

!$OMP master
                                         ! call fapp_stop("literm_other",1340,1)
                                           call clock_end(1340)
!$OMP end master


  END SUBROUTINE calc_dh

END MODULE GKV_advnc
