MODULE GKV_intgrl
!-------------------------------------------------------------------------------
!
!    Flux surface and field line averages
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public   intgrl_fsrf, intgrl_thet, intgrl_v0_moment, intgrl_v1_moment, intgrl_v2_moment, &
           intgrl_v0_moment_ms, intgrl_v2_moment_ms


    INTERFACE intgrl_fsrf
      module procedure  intgrl_fsrf_r, intgrl_fsrf_z
    END INTERFACE

    INTERFACE intgrl_thet
      module procedure  intgrl_thet_r, intgrl_thet_z
    END INTERFACE


CONTAINS


!note(fsrf does not require my dependence.)
!--------------------------------------
  SUBROUTINE intgrl_fsrf_r ( wn, wa )
!--------------------------------------
!     Flux-surface average of a real variable wn

    real(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: wn
    real(kind=DP), intent(out), &
      dimension(-nx:nx)                 :: wa

    real(kind=DP), dimension(-nx:nx)    :: ww
    integer  ::  mx, my, iz

          do mx = -nx, nx
            ww(mx)   = ( 0._DP, 0._DP )
          end do

      if ( rankw == 0 ) then


        my = 0
!$OMP parallel do reduction(+:ww)
          do iz = -nz, nz-1
            do mx = 0, nx
!              ww(mx)   = ww(mx) + wn(mx,my,iz) / omg(iz)
              ww(mx)   = ww(mx) + wn(mx,my,iz) * rootg(iz)
            end do
          end do


        call MPI_Allreduce( ww, wa, 2*nx+1, MPI_DOUBLE_PRECISION, &
                            MPI_SUM, zsp_comm_world, ierr_mpi )

!note(This procedure can be merged into upper one, 
!     since cfsrf does not have process dependence.)
        do mx = 0, nx
          wa(mx)   = wa(mx) / cfsrf
        end do

! --- reality condition
        do mx = 0, nx
          wa(-mx) = wa(mx)
        end do

      end if


  END SUBROUTINE intgrl_fsrf_r


!--------------------------------------
  SUBROUTINE intgrl_fsrf_z ( wn, wa )
!--------------------------------------
!     Flux-surface average of a complex variable wn

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1) :: wn
    complex(kind=DP), intent(out), &
      dimension(-nx:nx)                 :: wa

    complex(kind=DP), dimension(-nx:nx) :: ww
    integer  ::  mx, my, iz


          do mx = -nx, nx
            ww(mx)   = ( 0._DP, 0._DP )
          end do

      if ( rankw == 0 ) then

        my = 0
!$OMP parallel do reduction(+:ww)
          do iz = -nz, nz-1
            do mx = 0, nx
!              ww(mx)   = ww(mx) + wn(mx,my,iz) / omg(iz)
              ww(mx)   = ww(mx) + wn(mx,my,iz) * rootg(iz)
            end do
          end do


        call MPI_Allreduce( ww, wa, 2*nx+1, MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, zsp_comm_world, ierr_mpi )

        do mx = 0, nx
          wa(mx)   = wa(mx) / cfsrf
        end do

! --- reality condition
        do mx = 0, nx
          wa(-mx) = conjg( wa(mx) )
        end do

      end if


  END SUBROUTINE intgrl_fsrf_z


!--------------------------------------
  SUBROUTINE intgrl_thet_r ( wn, wa )
!--------------------------------------
!     average of a real variable wn in the theta space

    real(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1)     :: wn
    real(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)              :: wa

    real(kind=DP), dimension(-nx:nx,0:ny) :: ww
    integer  ::  mx, my, iz


!$OMP parallel do
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my)   = ( 0._DP, 0._DP )
          end do
        end do

!$OMP parallel do reduction(+:ww)
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
!            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) / omg(iz)
            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) * rootg(iz)
          end do
        end do
      end do

      call MPI_Allreduce( ww, wa, nxy, MPI_DOUBLE_PRECISION, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

      do my = ist_y, iend_y
        do mx = -nx, nx
          wa(mx,my)   = wa(mx,my) / cfsrf
        end do
      end do

! --- reality condition
      if ( rankw == 0 ) then
        my = 0
          do mx = 0, nx
            wa(-mx,my) = wa(mx,my)
          end do
      end if


  END SUBROUTINE intgrl_thet_r


!--------------------------------------
  SUBROUTINE intgrl_thet_z ( wn, wa )
!--------------------------------------
!     average of a complex variable wn in the theta space

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1)        :: wn
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny)                 :: wa

    complex(kind=DP), dimension(-nx:nx,0:ny) :: ww
    integer  ::  mx, my, iz


!$OMP parallel do
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my)   = ( 0._DP, 0._DP )
          end do
        end do

!$OMP parallel do reduction(+:ww)
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
!            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) / omg(iz)
            ww(mx,my)   = ww(mx,my) + wn(mx,my,iz) * rootg(iz)
          end do
        end do
      end do

      call MPI_Allreduce( ww, wa, nxy, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, zsp_comm_world, ierr_mpi )

!$OMP parallel do
      do my = ist_y, iend_y
        do mx = -nx, nx
          wa(mx,my)   = wa(mx,my) / cfsrf
        end do
      end do

! --- reality condition
      if ( rankw == 0 ) then
        my = 0
          do mx = 0, nx
            wa(-mx,my) = conjg( wa(mx,my) )
          end do
      end if


  END SUBROUTINE intgrl_thet_z


!--------------------------------------
  SUBROUTINE intgrl_v0_moment ( wf, wn )
!--------------------------------------
!     Calculate the zeroth order velocity moment of wf

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: wn

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ww
    real(kind=DP) :: cef
    integer :: mx, my, iz, iv, im

!$OMP parallel default(none) &
!$OMP shared(ww,wf,vp,dv,dvp,rankm,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,cef)

!$OMP do schedule(dynamic) 
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my,iz)   = ( 0._DP, 0._DP )
          end do
        end do
      end do
!$OMP end do 

      cef   = dv * twopi

      if( rankm == 0 ) then

        do im = 1, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                         &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

! for edge compensation
        im   = 1
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                           &
                        - ( - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) / 12._DP &
                          + ( wf(mx,my,iz,iv,im+1) * vp(iz,im+1)          &
                            - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) * 2._DP  &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do

      else

        do im = 0, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                         &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

      end if
!$OMP end parallel

      call MPI_Allreduce( ww, wn, nxyz, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, vel_comm_world, ierr_mpi )


  END SUBROUTINE intgrl_v0_moment


!--------------------------------------
  SUBROUTINE intgrl_v2_moment ( wf, wn )
!--------------------------------------
!     Calculate the second order velocity moment of wf

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: wn

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ww
    real(kind=DP) :: cef
    real(kind=DP) :: v2a, v2b
    integer :: mx, my, iz, iv, im

!$OMP parallel default(none) &
!$OMP shared(ww,wf,vl,vp,dv,dvp,rankm,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,cef,v2a,v2b)

!$OMP do schedule(dynamic) 
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my,iz)   = ( 0._DP, 0._DP )
          end do
        end do
      end do
!$OMP end do 

      cef   = dv * twopi

      if( rankm == 0 ) then

        do im = 1, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * ( vl(iv)**2 + vp(iz,im)**2 ) * 0.5_DP
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

! for edge compensation
        im   = 1
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  v2a   = ( vl(iv)**2 + vp(iz,im  )**2 ) * 0.5_DP
                  v2b   = ( vl(iv)**2 + vp(iz,im+1)**2 ) * 0.5_DP
                  ww(mx,my,iz)   = ww(mx,my,iz)                           &
                        - ( - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) / 12._DP &
                                * v2a                                     &
                          + ( wf(mx,my,iz,iv,im+1) * vp(iz,im+1)          &
                                * v2b                                     &
                            - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) * 2._DP  &
                                * v2a                                     &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do

      else

        do im = 0, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * ( vl(iv)**2 + vp(iz,im)**2 ) * 0.5_DP
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

      end if

!$OMP end parallel 

      call MPI_Allreduce( ww, wn, nxyz, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, vel_comm_world, ierr_mpi )

  END SUBROUTINE intgrl_v2_moment


!--------------------------------------
  SUBROUTINE intgrl_v1_moment ( wf, wn )
!--------------------------------------
!     Calculate the second order velocity moment of wf

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: wn

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ww
    real(kind=DP) :: cef
    integer :: mx, my, iz, iv, im

!$OMP parallel default(none) &
!$OMP shared(ww,wf,vl,vp,dv,dvp,rankm,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,cef)

!$OMP do schedule(dynamic) 
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my,iz)   = ( 0._DP, 0._DP )
          end do
        end do
      end do
!$OMP end do 

      cef   = dv * twopi

      if( rankm == 0 ) then

        do im = 1, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * vl(iv) 
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

! for edge compensation
        im   = 1
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                           &
                        - ( - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) / 12._DP &
                                * vl(iv)                                     &
                          + ( wf(mx,my,iz,iv,im+1) * vp(iz,im+1)          &
                                * vl(iv)                                     &
                            - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) * 2._DP  &
                                * vl(iv)                                     &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do

      else

        do im = 0, nm
          do iv = 1, 2*nv
!%%%% do schedule(dynamic) 
!$OMP do
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * vl(iv)
                end do
              end do
            end do
!%%%% end do nowait
!$OMP end do
          end do
        end do

      end if

!$OMP end parallel

      call MPI_Allreduce( ww, wn, nxyz, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, vel_comm_world, ierr_mpi )

  END SUBROUTINE intgrl_v1_moment


!--------------------------------------
  SUBROUTINE intgrl_v0_moment_ms ( wf, wn )
!--------------------------------------
!     Calculate the zeroth order velocity moment of wf over species

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: wn

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ww
    real(kind=DP) :: cef
    integer :: mx, my, iz, iv, im

!$OMP parallel default(none) &
!$OMP shared(ww,wf,vp,dv,dvp,rankm,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,cef)

!$OMP do schedule(dynamic) 
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my,iz)   = ( 0._DP, 0._DP )
          end do
        end do
      end do
!$OMP end do 

      cef   = dv * twopi

      if( rankm == 0 ) then

        do im = 1, nm
!$OMP do schedule(dynamic) reduction(+:ww)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                         &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do
          end do
!$OMP end do 
        end do

! for edge compensation
        im   = 1
!$OMP do schedule(dynamic) reduction(+:ww)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                           &
                        - ( - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) / 12._DP &
                          + ( wf(mx,my,iz,iv,im+1) * vp(iz,im+1)          &
                            - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) * 2._DP  &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)
                end do
              end do
            end do
          end do
!$OMP end do 

      else

!$OMP do schedule(dynamic) reduction(+:ww)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                         &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef
                end do
              end do
            end do
          end do
        end do
!$OMP end do 

      end if

!$OMP end parallel

      call MPI_Allreduce( ww, wn, nxyz, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, spc_comm_world, ierr_mpi )



  END SUBROUTINE intgrl_v0_moment_ms


!--------------------------------------
  SUBROUTINE intgrl_v2_moment_ms ( wf, wn )
!--------------------------------------
!     Calculate the second order velocity moment of wf over species

    complex(kind=DP), intent(in),  &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)   :: wf
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: wn

    complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1) :: ww
    real(kind=DP) :: cef
    real(kind=DP) :: v2a, v2b
    integer :: mx, my, iz, iv, im

!$OMP parallel default(none) &
!$OMP shared(ww,wf,vl,vp,dv,dvp,rankm,ist_y,iend_y) &
!$OMP private(mx,my,iz,iv,im,cef,v2a,v2b)

!$OMP do schedule(dynamic) 
      do iz = -nz, nz-1
        do my = ist_y, iend_y
          do mx = -nx, nx
            ww(mx,my,iz)   = ( 0._DP, 0._DP )
          end do
        end do
      end do
!$OMP end do 

      cef   = dv * twopi

      if( rankm == 0 ) then

!$OMP do schedule(dynamic) reduction(+:ww) 
        do im = 1, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * ( vl(iv)**2 + vp(iz,im)**2 ) * 0.5_DP
                end do
              end do
            end do
          end do
        end do
!$OMP end do 

! for edge compensation
        im   = 1
!$OMP do schedule(dynamic) reduction(+:ww) 
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  v2a   = ( vl(iv)**2 + vp(iz,im  )**2 ) * 0.5_DP
                  v2b   = ( vl(iv)**2 + vp(iz,im+1)**2 ) * 0.5_DP
                  ww(mx,my,iz)   = ww(mx,my,iz)                           &
                        - ( - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) / 12._DP &
                                * v2a                                     &
                          + ( wf(mx,my,iz,iv,im+1) * vp(iz,im+1)          &
                                * v2b                                     &
                            - wf(mx,my,iz,iv,im  ) * vp(iz,im  ) * 2._DP  &
                                * v2a                                     &
                            ) * 11._DP / 720._DP                          &
                          ) * cef * dvp(iz)
                end do
              end do
            end do
          end do
!$OMP end do 

      else

!$OMP do schedule(dynamic) reduction(+:ww) 
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  ww(mx,my,iz)   = ww(mx,my,iz)                          &
                        + wf(mx,my,iz,iv,im) * vp(iz,im) * dvp(iz) * cef &
                          * ( vl(iv)**2 + vp(iz,im)**2 ) * 0.5_DP
                end do
              end do
            end do
          end do
        end do
!$OMP end do 

      end if

!$OMP end parallel

      call MPI_Allreduce( ww, wn, nxyz, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, spc_comm_world, ierr_mpi )


  END SUBROUTINE intgrl_v2_moment_ms


END MODULE GKV_intgrl
