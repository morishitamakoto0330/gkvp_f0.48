MODULE GKV_dtc
!-------------------------------------------------------------------------------
!
!    Time step size control
!
!      GKV-plus f0.25 ( S.Maeyama, Nov 2012)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_exb, only: exw_xw, eyw_xw, bxw_xw, byw_xw

  implicit none

  private

  real(kind=DP), save :: dt_linear, dt_nl, dt_limit

  real(kind=DP), save :: dx_inv, dy_inv

  public   dtc_init, dtc_cntrl


CONTAINS


!--------------------------------------
  SUBROUTINE dtc_init( lx, ly, vmax )
!--------------------------------------

    real(kind=DP), intent(in) :: lx, ly, vmax

    real(kind=DP) :: dt_perp, dt_zz, dt_vl, dt_col
    real(kind=DP) :: kvd_max, kvd_max2, vl_max, vl_max2, mir_max, mir_max2
    real(kind=DP) :: ksq_max0, ksq_max, nu_max, nu_max2, nu_temp
    real(kind=DP) :: cs, dx, dy
    integer :: mx, my, iz, iv, im


      ksq_max0 = 0._DP
      ksq_max  = 0._DP
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
!              do mx = 0, nx
                if ( ksq_max0 < ksq(mx,my,iz) ) ksq_max0 = ksq(mx,my,iz)
              end do
            end do
          end do
      call MPI_Allreduce( ksq_max0, ksq_max, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )

      kvd_max = 0._DP
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
!              do mx = 0, nx
                if ( kvd_max < kvd(mx,my,iz,iv,im) ) then 
                  kvd_max = kvd(mx,my,iz,iv,im)
                end if 
              end do
            end do
          end do
        end do
      end do
      call MPI_Allreduce( kvd_max, kvd_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      dt_perp = courant_num * pi / kvd_max2

      cs = sqrt( tau(ranks) / Anum(ranks) )
      vl_max = 0._DP
      do iz = -nz, nz-1
        if ( vl_max < cs * vmax / dpara(iz) ) vl_max = cs * vmax / dpara(iz)
      end do
      call MPI_Allreduce( vl_max, vl_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      dt_zz = courant_num / vl_max2

      mir_max = 0._DP
      do im = 0, nm
        do iz = -nz, nz-1
          if ( mir_max < cs * mir(iz,im) ) mir_max = cs * mir(iz,im)
        end do
      end do
      call MPI_Allreduce( mir_max, mir_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
      dt_vl = courant_num * dv / mir_max2

      if ( trim(col_type) == "LB" ) then

        nu_max = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP &
               * ( 2._DP / dv**2 )
        do iz = -nz, nz-1
          nu_temp = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP &
                  * ( 2._DP / dvp(iz)**2 )
          if ( nu_max < nu_temp ) nu_max = nu_temp
        end do
        call MPI_Allreduce( nu_max, nu_max2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
        dt_col = courant_num / nu_max2
        
      else if ( trim(col_type) == "full" ) then

        nu_max = 0._DP
        do im = 0, nm
          do iv = 1, nv
            do iz = -nz, nz-1
              nu_temp = ( nu_ps(iz,iv,im) * vl(iv)**2    & 
                        + nu_ds(iz,iv,im) * vp(iz,im)**2 &
                        ) * 0.5_DP * ( 2._DP / dv**2 )
              if ( nu_max < nu_temp ) nu_max = nu_temp
            end do
          end do
        end do
        do im = 0, nm
          do iv = 1, nv
            do iz = -nz, nz-1
              nu_temp = ( nu_ds(iz,iv,im) * vl(iv)**2    & 
                        + nu_ps(iz,iv,im) * vp(iz,im)**2 &
                        ) * 0.5_DP * ( 2._DP / dvp(iz)**2 )
              if ( nu_max < nu_temp ) nu_max = nu_temp
            end do
          end do
        end do
        call MPI_Allreduce( nu_max, nu_max2, 1, MPI_DOUBLE_PRECISION, &
                            MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
        dt_col = courant_num / nu_max2

      else

        write( olog, * ) "This col_type is not supported by dtc:", trim(col_type)
        dt_col = 99._DP

      end if


      dt_linear = min( dt_perp, dt_zz, dt_vl, dt_col )

      dt_limit = min( dt_max, dt_linear )

      dt = dt_max

      if ( adapt_dt ) then
        if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) dt = dt_limit
      end if

        write( olog, * ) " # Time step size control"
        write( olog, * ) ""
        write( olog, * ) " # courant num. = ", courant_num
        write( olog, * ) " # ksq_max      = ", ksq_max
        write( olog, * ) " # dt_perp      = ", dt_perp
        write( olog, * ) " # dt_zz        = ", dt_zz
        write( olog, * ) " # dt_vl        = ", dt_vl
        write( olog, * ) " # dt_col       = ", dt_col
        write( olog, * ) " # dt_linear    = ", dt_linear
        write( olog, * ) " # dt_max       = ", dt_max
        write( olog, * ) " # dt           = ", dt
        write( olog, * ) ""

      dx = lx / real( nxw, kind=DP )
      dy = ly / real( nyw, kind=DP )
      dx_inv = 1._DP / dx
      dy_inv = 1._DP / dy

  
  END SUBROUTINE dtc_init


!--------------------------------------
  SUBROUTINE dtc_cntrl( time, id )
!--------------------------------------

    real(kind=DP), intent(in) :: time
    integer, intent(in) :: id

    real(kind=DP), save :: tout_dtc


      if( id == 0 ) then

        tout_dtc  = ( int( ( time + eps )/dtout_dtc ) + 1 ) * dtout_dtc
 
        call dtc_estimate

        if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) then
          dt = dt_limit
          write( olog, * ) &
            " # dt is changed at time = ", time, ", dt = ", dt
        end if
     
        if ( rankg == 0 ) then
          write( unit=odtc, fmt="(f10.5, 1p, 3e15.7)" )  &
            time, dt, dt_limit, dt_nl
        end if
  
      else if( id == 1 ) then

        if ( time >= tout_dtc - eps ) then

          tout_dtc   = tout_dtc + dtout_dtc
   
          call dtc_estimate
  
          if ( dt < 0.8_DP * dt_limit .or. 1.2_DP * dt_limit < dt ) then
            dt = dt_limit
            write( olog, * ) &
              " # dt is changed at time = ", time, ", dt = ", dt
          end if
       
          if ( rankg == 0 ) then
            write( unit=odtc, fmt="(f10.5, 1p, 3e15.7)" )  &
              time, dt, dt_limit, dt_nl
          end if

        end if

      end if


  END SUBROUTINE dtc_cntrl


!--------------------------------------
  SUBROUTINE dtc_estimate
!--------------------------------------

    real(kind=DP) :: wx_nl1, wx_nl2, wy_nl1, wy_nl2, w_nl_max, w_nl_max2
    real(kind=DP) :: cs
    integer  ::  mx, my, iz, iv, im
  
  
      w_nl_max = 0._DP
  
      cs = sqrt( tau(ranks) / Anum(ranks) )

      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do mx = ist_xw, iend_xw
              do my = 0, nyw
                wx_nl1 = abs( dx_inv * real( eyw_xw(my,mx,iz,im) - cs * vl(iv) &
                                             * byw_xw(my,mx,iz,im), kind=DP ) )
                wx_nl2 = abs( dx_inv * aimag( eyw_xw(my,mx,iz,im) - cs * vl(iv) &
                                             * byw_xw(my,mx,iz,im)          ) )
                wy_nl1 = abs( dy_inv * real( exw_xw(my,mx,iz,im) - cs * vl(iv) &
                                             * bxw_xw(my,mx,iz,im), kind=DP ) )
                wy_nl2 = abs( dy_inv * aimag( exw_xw(my,mx,iz,im) - cs * vl(iv) &
                                             * bxw_xw(my,mx,iz,im)          ) )
                if ( w_nl_max < wx_nl1 ) w_nl_max = wx_nl1
                if ( w_nl_max < wx_nl2 ) w_nl_max = wx_nl2
                if ( w_nl_max < wy_nl1 ) w_nl_max = wy_nl1
                if ( w_nl_max < wy_nl2 ) w_nl_max = wy_nl2
              end do
            end do
          end do
        end do
      end do
  
      call MPI_Allreduce( w_nl_max, w_nl_max2, 1, MPI_DOUBLE_PRECISION, &
                          MPI_MAX, MPI_COMM_WORLD, ierr_mpi )
  
      dt_nl = courant_num / w_nl_max2
  
      dt_limit = min( dt_max, dt_linear, dt_nl )

  
  END SUBROUTINE dtc_estimate


END MODULE GKV_dtc
