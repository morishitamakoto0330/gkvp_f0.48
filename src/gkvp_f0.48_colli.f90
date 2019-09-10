MODULE GKV_colli
!-------------------------------------------------------------------------------
!
!    Collision term
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_clock, only : clock_sta, clock_end

  implicit none

  private

  public   colli_LB_model, colli_GK_CT, colli_GK_CT6, colli_GK_CF_DT, & 
           colli_set_param, colli_moment_calc, colli_moment_redc, colli_comm_alltoall, & 
           colli_dfdvp, colli_dfdvp6, colli_zeroset, colli_hhset, colli_wwset


CONTAINS


!--------------------------------------
  SUBROUTINE colli_set_param (q0, eps_r, nust)
!-------------------------------------------------------------------------------
!
!    Set parameters for GK collision term
!
!    by M. Nakata and M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    real(kind=DP), parameter :: mp      = 1.67262178d-24, & ! proton mass in g
                                ee      = 4.80320425d-10, & ! elementary charge in esu
                                ev2erg  = 1.60217657d-12    ! erg/eV  (cf. 1J = 10^7 erg)
    
    real(kind=DP),                    intent(in)  :: q0, eps_r
    real(kind=DP), dimension(0:ns-1,0:ns-1), intent(out) :: nust

    real(kind=DP), dimension(0:ns-1)        :: tmpr, dens, freq_factor
    real(kind=DP), dimension(0:ns-1,0:ns-1) :: log_lambda, calpha, ctheta, cgamma, ceta, cxi
    real(kind=DP)                           :: cph, dph, cgg

    integer :: is, is1, is2, iz, iv, im, mx, my
  

! --- temperature [in eV] and density [in cm^(-3)]
    do is = 0, ns-1
      tmpr(is) = tau(is) * Tref*1.d3
      dens(is) = Nref*1.d-6 * fcs(is)/Znum(is)
    enddo

! --- factor for collision frequencies
    do is = 0, ns-1
    freq_factor(is)  = (dens(is) * ee**4 * Lref*1.d2) / (Tref*1.d3*ev2erg)**2 
    end do


! --- Coulomb logarithm in cm^(-3) and eV units (see NRL plasma Formulary)  
    do is1 = 0, ns-1
      if (sgn(is1) < 0.d0) then  !! For is1 = electron

        do is2 = 0, ns-1 
          if (sgn(is2) < 0.d0) then   !! e-e case
            log_lambda(is1,is2) = 23.5_DP - dlog( dsqrt( dens(is1) ) * tmpr(is1)**(-1.25_DP) )  &
                                          - dsqrt( 1.d-5 + (( dlog(tmpr(is1)) - 2._DP )**2 )/16._DP )

          else                        !! e-i case
            log_lambda(is1,is2) = 24._DP - dlog( dsqrt( dens(is1) ) / tmpr(is1) )
          endif
        enddo

      else                     !! For is1 = ions

        do is2 = 0, ns-1 
          if (sgn(is2) < 0.d0) then   !! i-e case
            log_lambda(is1,is2) = 24._DP - dlog( dsqrt( dens(is2) ) / tmpr(is2) )

          else                       !! i-i case
            log_lambda(is1,is2) = &
              23._DP - dlog( Znum(is1)*Znum(is2)*(Anum(is1)+Anum(is2))/(Anum(is1)*tmpr(is2)+Anum(is2)*tmpr(is1)) &
                             * dsqrt( (dens(is1) * Znum(is1)**2)/tmpr(is1)                                       &
                                     + (dens(is2) * Znum(is2)**2)/tmpr(is2) ) )
          endif
        enddo

      endif
    enddo

! --- Constant parameters
    do is1 = 0, ns-1
      do is2 = 0, ns-1 

         ctauiv(is1,is2) = freq_factor(is2) * (8._DP*dsqrt(pi)/3._DP/dsqrt(2._DP))*log_lambda(is1,is2)  & 
                                   * ( Znum(is1)**2*Znum(is2)**2/dsqrt(Anum(is1))/tau(is1)**1.5 )

         calpha(is1,is2) = dsqrt( tau(is1) * Anum(is2) / ( tau(is2) * Anum(is1) ) )
         ctheta(is1,is2) = dsqrt( tau(is1) * ( Anum(is1) + Anum(is2) ) / ( tau(is1) * Anum(is2) + tau(is2) * Anum(is1) ) )

         cgamma(is1,is2) = - Anum(is1) * calpha(is1,is2)                                      &
                            * (tau(is1)/tau(is2) + calpha(is1,is2)**2) * ctauiv(is1,is2)      &
                             / (1._DP + calpha(is1,is2)**2)**1.5_DP 

         ceta(is1,is2)   = - tau(is1) * 3._DP * calpha(is1,is2)                               &
                            * (tau(is1)/tau(is2) + calpha(is1,is2)**2) * ctauiv(is1,is2)      &
                             / (1._DP + calpha(is1,is2)**2)**2.5_DP

          cxi(is1,is2)   =  calpha(is1,is2) * ( ctheta(is1,is2) - 1._DP ) * ctauiv(is1,is2)   &
                             / dsqrt(1._DP + calpha(is1,is2)**2) 

         nust(is1,is2)   = q0*(ctauiv(is1,is2)/dsqrt(2._DP))/(eps_r**1.5*dsqrt(tau(is1)/Anum(is1)))

      enddo
    enddo

! --- xxa = v/vta/sqrt(2), where vta = sqrt(Ta/ma)
    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1
          xxa(iz,iv,im) = dsqrt(vl(iv)**2 + vp(iz,im)**2)/dsqrt(2._DP) 
        end do 
      end do
    end do

! --- collision frequencies and v-space functions
    do im = 0, nm
      do iv = 1, 2*nv
        do iz = -nz, nz-1 
          do is1 = 0, ns-1
            do is2 = 0, ns-1 

              cph = derf(calpha(is1,is2)*xxa(iz,iv,im))
              dph = 2._DP/dsqrt(pi)*dexp(-calpha(is1,is2)**2*xxa(iz,iv,im)**2)
              cgg = (cph - calpha(is1,is2)*xxa(iz,iv,im)*dph)/(calpha(is1,is2)**2*xxa(iz,iv,im)**2)*0.5_DP

              nu_d(iz,iv,im,is1,is2) = 0.75_DP*dsqrt(pi)*ctauiv(is1,is2)*(cph-cgg)/xxa(iz,iv,im)**3
              nu_p(iz,iv,im,is1,is2) = 1.50_DP*dsqrt(pi)*ctauiv(is1,is2)*(  cgg  )/xxa(iz,iv,im)**3
              nu_h(iz,iv,im,is1,is2) = 0.75_DP*dsqrt(pi)*ctauiv(is1,is2)*calpha(is1,is2)*dph/xxa(iz,iv,im)**2
              nu_g(iz,iv,im,is1,is2) = nu_p(iz,iv,im,is1,is2)*xxa(iz,iv,im)**2*(1._DP-calpha(is1,is2)**2)

              c_t0(iz,iv,im,is1,is2,1)  = - (1._DP + calpha(is1,is2)**2)*fmx(iz,iv,im)*nu_p(iz,iv,im,is1,is2)               &
                                             * xxa(iz,iv,im)**2*vl(iv)
              c_t0(iz,iv,im,is1,is2,2)  = - 1.5_DP*dsqrt(pi)*ctauiv(is1,is2)*fmx(iz,iv,im)                                  & 
                                             * ( cph - calpha(is1,is2)*xxa(iz,iv,im)*(1._DP + calpha(is1,is2)**2)*dph )     & 
                                             / calpha(is1,is2)**2 / xxa(iz,iv,im)

              x_tst(1,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*fmx(iz,iv,im)*vl(iv)
              x_tst(2,iz,iv,im,is1,is2) = x_tst(1,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv) 
              x_tst(3,iz,iv,im,is1,is2) = x_tst(1,iz,iv,im,is1,is2)*(xxa(iz,iv,im)**2/1.5_DP - 1._DP)/vl(iv)
              x_tst(4,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*( c_t0(iz,iv,im,is1,is2,1)                              &
                                                      - (ctheta(is1,is2) - 1._DP)*calpha(is1,is2)*ctauiv(is1,is2)           &
                                                              * fmx(iz,iv,im)*vl(iv)/dsqrt(1._DP + calpha(is1,is2)**2) )
              x_tst(5,iz,iv,im,is1,is2) =  x_tst(4,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv)  
              x_tst(6,iz,iv,im,is1,is2) = (ctheta(is1,is2) - 1._DP)*( c_t0(iz,iv,im,is1,is2,2)*2._DP/3._DP                &
                                                      - (ctheta(is1,is2) - 1._DP)*calpha(is1,is2)*ctauiv(is1,is2)         &
                                                              * fmx(iz,iv,im)*(xxa(iz,iv,im)**2/1.5_DP - 1._DP)*2._DP     &
                                                               / (1._DP + calpha(is1,is2)**2)**1.5 )

              y_fld(1,iz,iv,im,is1,is2) = - (fcs(is2)/Znum(is2))/(fcs(is1)/Znum(is1))*calpha(is1,is2)*Anum(is1)             & 
                                                    * tau(is2)*ctheta(is1,is2)*ctheta(is2,is1)/tau(is1)/cgamma(is1,is2)     &
                                                    * ( c_t0(iz,iv,im,is1,is2,1) - cxi(is1,is2)*fmx(iz,iv,im)*vl(iv) ) 
              y_fld(2,iz,iv,im,is1,is2) = y_fld(1,iz,iv,im,is1,is2)*vp(iz,im)/vl(iv)  
              y_fld(3,iz,iv,im,is1,is2) = - (fcs(is2)/Znum(is2))/(fcs(is1)/Znum(is1))                                     & 
                                                    * tau(is2)*ctheta(is1,is2)*ctheta(is2,is1)/ceta(is1,is2)              &
                                                    * ( c_t0(iz,iv,im,is1,is2,2)                                          &
                                                           - cxi(is1,is2)/(1._DP+calpha(is1,is2)**2)*fmx(iz,iv,im)        &
                                                               *(2._DP*xxa(iz,iv,im)**2 - 3._DP) ) 
              y_fld(4,iz,iv,im,is1,is2) = - y_fld(1,iz,iv,im,is1,is2)*cxi(is2,is1) 
              y_fld(5,iz,iv,im,is1,is2) = - y_fld(2,iz,iv,im,is1,is2)*cxi(is2,is1) 
              y_fld(6,iz,iv,im,is1,is2) = - y_fld(3,iz,iv,im,is1,is2)*2._DP*cxi(is2,is1)/(1._DP+calpha(is2,is1)**2) 

            end do 
          end do 
        end do 
      end do
    end do 

! --- summation of collision frequencies with respect to is2, and adiabatic term (used in colli_GK_CT)
    nu_hs = 0._DP 
    nu_gs = 0._DP 
    nu_ds = 0._DP 
    nu_ps = 0._DP 
    is1 = ranks
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1 
            do is2 = 0, ns-1 
                nu_hs(iz,iv,im) = nu_hs(iz,iv,im) + nu_h(iz,iv,im,is1,is2)
                nu_gs(iz,iv,im) = nu_gs(iz,iv,im) + nu_g(iz,iv,im,is1,is2)
                nu_ds(iz,iv,im) = nu_ds(iz,iv,im) + nu_d(iz,iv,im,is1,is2)
                nu_ps(iz,iv,im) = nu_ps(iz,iv,im) + nu_p(iz,iv,im,is1,is2)
            end do 
          end do 
        end do
      end do 


! --- adiabatic part (used in colli_GK_CT)
    is1 = ranks
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1 
            do my = ist_y, iend_y 
              do mx = 0, nx
                adbtc(mx,my,iz,iv,im) =                                              &
                               ( -( nu_ds(iz,iv,im)*vl(iv)**2                        & 
                                     + nu_ps(iz,iv,im)*vp(iz,im)**2 )                & 
                                   * 0.25_DP*ksq(mx,my,iz)*Anum(is1)                 &
                                           /Znum(is1)/omg(iz)**2                     &
                                   * (j0(mx,my,iz,im) - j2(mx,my,iz,im))             &
                                 -( nu_hs(iz,iv,im)*vp(iz,im)                        &
                                     - 0.5_DP*nu_ps(iz,iv,im)*vp(iz,im)              &
                                             *(1._DP-vl(iv)**2-vp(iz,im)**2)         &
                                     + 0.5_DP*nu_ds(iz,iv,im)                        &
                                             *(vl(iv)**2/vp(iz,im)-vp(iz,im)) )      &
                                   * dsqrt(ksq(mx,my,iz)*Anum(is1)/tau(is1))/omg(iz) & 
                                   * j1(mx,my,iz,im)                                 &
                                 -( nu_ds(iz,iv,im)                                  &
                                             *(2._DP*vl(iv)**2+vp(iz,im)**2)         & 
                                     + nu_ps(iz,iv,im)*vp(iz,im)**2 )                & 
                                   * 0.25_DP*ksq(mx,my,iz)*Anum(is1)                 &
                                           /Znum(is1)/omg(iz)**2 * j0(mx,my,iz,im)   &          
                               ) * fmx(iz,iv,im)*sgn(is1)*real(iFLR, kind=DP)

                adbtc0(mx,my,iz,iv,im) =                                             &
                               ( nu_ds(iz,iv,im)*vl(iv)**2                           & 
                                   * ksq(mx,my,iz)*Anum(is1)/Znum(is1)/omg(iz)**2    &
                               ) * fmx(iz,iv,im)*sgn(is1)*real(iFLR, kind=DP)
              end do 
            end do 
          end do 
        end do 
      end do


! --- set v-space functions used in colli_moment

    if ( iFLR == 1 ) then
      is1 = ranks
        do is2 = 0, ns-1
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx,nx
                    vfunc(mx,my,iz,iv,im,is2,1) = & 
                                    j0(mx,my,iz,im)*c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
                    vfunc(mx,my,iz,iv,im,is2,2) = & 
                                    j1(mx,my,iz,im)*vp(iz,im)*c_t0(iz,iv,im,is1,is2,1)/vl(iv)/fmx(iz,iv,im)
                    vfunc(mx,my,iz,iv,im,is2,3) = & 
                                    j0(mx,my,iz,im)*c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
                    vfunc(mx,my,iz,iv,im,is2,4) = j0(mx,my,iz,im)*vl(iv)
                    vfunc(mx,my,iz,iv,im,is2,5) = j1(mx,my,iz,im)*vp(iz,im)
                    vfunc(mx,my,iz,iv,im,is2,6) = j0(mx,my,iz,im)*(xxa(iz,iv,im)**2-1.5_DP)
                  end do 
                end do 
              end do 
            end do 
          end do 
        end do 
   
    else 

      is1 = ranks
        do is2 = 0, ns-1
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx,nx
                    vfunc(mx,my,iz,iv,im,is2,1) = c_t0(iz,iv,im,is1,is2,1)/fmx(iz,iv,im) 
                    vfunc(mx,my,iz,iv,im,is2,2) = 0._DP
                    vfunc(mx,my,iz,iv,im,is2,3) = c_t0(iz,iv,im,is1,is2,2)/fmx(iz,iv,im)
                    vfunc(mx,my,iz,iv,im,is2,4) = vl(iv)
                    vfunc(mx,my,iz,iv,im,is2,5) = 0._DP
                    vfunc(mx,my,iz,iv,im,is2,6) = (xxa(iz,iv,im)**2-1.5_DP)
                  end do 
                end do 
              end do 
            end do 
          end do 
        end do 

    end if


! -----------------------------------
! --- Output constants
    if ( rankg == nprocz/2 ) then

      do is1 = 0, ns-1
        do is2 = 0, ns-1 
          write(unit=ocst,fmt="(2I3,SP,256ES24.15e3)") is1, is2, ctheta(is1,is2), calpha(is1,is2), &
                                                                 fcs(is1)/Znum(is1)*ceta(is1,is2), &
                                                               fcs(is1)/Znum(is1)*cgamma(is1,is2), & 
                                                                    cxi(is1,is2), ctauiv(is1,is2), &
                                                                              log_lambda(is1,is2)
! --- Note that, for ns >=3, cgamma(is1,is2) /= cgamma(is2,is1), but dens(is1)*cgamma(is1,is2) = dense(is2)*cgamma(is2,is1)
! ---  due to normalizartion with dens(is). 
        enddo
      enddo

!! --- for debug
!      iz = -nz
!      do im = 0, nm
!        do iv = 1, 2*nv
!          write(unit=4000,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_h(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_g(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_d(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((nu_p(iz,iv,im,is1,is2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                  ((c_t0(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=4005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                  ((c_t0(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!        end do
!        write (unit=4000,fmt=*)
!        write (unit=4001,fmt=*)
!        write (unit=4002,fmt=*)
!        write (unit=4003,fmt=*)
!        write (unit=4004,fmt=*)
!        write (unit=4005,fmt=*)
!      end do
!
!! --- for debug
!      iz = -nz
!      do im = 0, nm
!        do iv = 1, 2*nv
!          write(unit=5001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,3), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,4), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,5), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=5006,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((x_tst(iz,iv,im,is1,is2,6), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6001,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,1), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6002,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,2), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6003,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,3), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6004,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,4), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6005,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,5), is2 = 0, ns-1), is1 = 0, ns-1)
!          write(unit=6006,fmt="(2f15.8,SP,256ES24.15e3)") vl(iv), vp(iz,im), & 
!                                    ((y_fld(iz,iv,im,is1,is2,6), is2 = 0, ns-1), is1 = 0, ns-1)
!        end do
!        write (unit=5001,fmt=*)
!        write (unit=5002,fmt=*)
!        write (unit=5003,fmt=*)
!        write (unit=5004,fmt=*)
!        write (unit=5005,fmt=*)
!        write (unit=5006,fmt=*)
!        write (unit=6001,fmt=*)
!        write (unit=6002,fmt=*)
!        write (unit=6003,fmt=*)
!        write (unit=6004,fmt=*)
!        write (unit=6005,fmt=*)
!        write (unit=6006,fmt=*)
!      end do

    end if

    return

   END SUBROUTINE colli_set_param


!--------------------------------------
  SUBROUTINE colli_LB_model( ff, im, cf )
!--------------------------------------
!   Lenard-Bernstein model collsion operator

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    integer, intent(in) :: im
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    real(kind=DP) :: nu_s, cef1, cef2
    real(kind=DP), dimension(-nz:nz-1) :: cef3, cef4
    integer  ::  mx, my, iz, iv


!$OMP master
                                           call clock_sta(1311)
                                         ! call fapp_start("literm_colli_ct",1311,1)
!$OMP end master

! --- Note that nu(ranks) is a bias factor 
      nu_s = nu(ranks)*3._DP*dsqrt(pi)*ctauiv(ranks,ranks)/4._DP
!      nu_s = 1.d-3

      cef1   = nu_s / ( 12._DP * dv * dv )
      cef2   = nu_s / ( 12._DP * dv )

      do iz = -nz, nz-1
        cef3(iz)   = nu_s / ( 12._DP * dvp(iz) * dvp(iz) )
        cef4(iz)   = nu_s / ( 12._DP * dvp(iz) )
      end do

      if( rankm /= 0  ) then

          do iv = 1, 2*nv
!$OMP do schedule (dynamic)
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  cf(mx,my,iz,iv) =                                           &
                            ( -          ff(mx,my,iz,iv+2,im)                 &
                              + 16._DP * ff(mx,my,iz,iv+1,im)                 &
                              - 30._DP * ff(mx,my,iz,iv  ,im)                 &
                              + 16._DP * ff(mx,my,iz,iv-1,im)                 &
                              -          ff(mx,my,iz,iv-2,im)                 &
                             ) * cef1                                         &
                           + ( -          ff(mx,my,iz,iv+2,im)                &
                               +  8._DP * ff(mx,my,iz,iv+1,im)                &
                               -  8._DP * ff(mx,my,iz,iv-1,im)                &
                               +          ff(mx,my,iz,iv-2,im)                &
                             ) * cef2 * vl(iv)                                &
                           + ( -          ff(mx,my,iz,iv,im+2)                &
                               + 16._DP * ff(mx,my,iz,iv,im+1)                &
                               - 30._DP * ff(mx,my,iz,iv,im  )                &
                               + 16._DP * ff(mx,my,iz,iv,im-1)                &
                               -          ff(mx,my,iz,iv,im-2)                &
                             ) * cef3(iz)                                     &
                           + ( -          ff(mx,my,iz,iv,im+2)                &
                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
                               +          ff(mx,my,iz,iv,im-2)                &
                             ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                           + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &      
                           - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
                             / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
                               * real(iFLR, kind=DP)
                end do
              end do
            end do
!$OMP end do nowait
          end do

      else

          if ( im == 0 ) then

            do iv = 1, 2*nv
!$OMP do schedule (dynamic)
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    cf(mx,my,iz,iv) =                                           &
                               ( -          ff(mx,my,iz,iv+2,im)                &
                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
                                 -          ff(mx,my,iz,iv-2,im)                &
                               ) * cef1                                         &
                             + ( -          ff(mx,my,iz,iv+2,im)                &
                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
                                 +          ff(mx,my,iz,iv-2,im)                &
                               ) * cef2 * vl(iv)                                &
                             + ( -          ff(mx,my,iz,iv,im+2)                &
                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
                                 -          ff(mx,my,iz,iv,im+2)                &
                               ) * cef3(iz) * 2._DP                             &
                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
                               * real(iFLR, kind=DP)
                  end do
                end do
              end do
!$OMP end do nowait
            end do

          else if ( im == 1 ) then

            do iv = 1, 2*nv
!$OMP do schedule (dynamic)
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    cf(mx,my,iz,iv) =                                           &
                               ( -          ff(mx,my,iz,iv+2,im)                &
                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
                                 -          ff(mx,my,iz,iv-2,im)                &
                               ) * cef1                                         &
                             + ( -          ff(mx,my,iz,iv+2,im)                &
                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
                                 +          ff(mx,my,iz,iv-2,im)                &
                               ) * cef2 * vl(iv)                                &
                             + ( -          ff(mx,my,iz,iv,im+2)                &
                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
                                 + 16._DP * ff(mx,my,iz,iv,im-1)                &
                                 -          ff(mx,my,iz,iv,im  )                &
                               ) * cef3(iz)                                     &
                             + ( -          ff(mx,my,iz,iv,im+2)                &
                                 +  8._DP * ff(mx,my,iz,iv,im+1)                &
                                 -  8._DP * ff(mx,my,iz,iv,im-1)                &
                                 +          ff(mx,my,iz,iv,im  )                &
                               ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &   
                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) &
                               * real(iFLR, kind=DP)
                  end do
                end do
              end do
!$OMP end do nowait
            end do

          else  ! 2=<im=<nm

            do iv = 1, 2*nv
!$OMP do schedule (dynamic)
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    cf(mx,my,iz,iv) =                                           &
                               ( -          ff(mx,my,iz,iv+2,im)                &
                                 + 16._DP * ff(mx,my,iz,iv+1,im)                &
                                 - 30._DP * ff(mx,my,iz,iv  ,im)                &
                                 + 16._DP * ff(mx,my,iz,iv-1,im)                &
                                 -          ff(mx,my,iz,iv-2,im)                &
                               ) * cef1                                         &
                             + ( -          ff(mx,my,iz,iv+2,im)                &
                                 +  8._DP * ff(mx,my,iz,iv+1,im)                &
                                 -  8._DP * ff(mx,my,iz,iv-1,im)                &
                                 +          ff(mx,my,iz,iv-2,im)                &
                               ) * cef2 * vl(iv)                                &
                             + ( -          ff(mx,my,iz,iv,im+2)                &
                                 + 16._DP * ff(mx,my,iz,iv,im+1)                &
                                 - 30._DP * ff(mx,my,iz,iv,im  )                &
                                 + 16._DP * ff(mx,my,iz,iv,im-1)                &
                                 -          ff(mx,my,iz,iv,im-2)                &
                               ) * cef3(iz)                                     &
                             + ( -          ff(mx,my,iz,iv,im+2)                &
                                 +  8._DP * ff(mx,my,iz,iv,im+1)                &
                                 -  8._DP * ff(mx,my,iz,iv,im-1)                &
                                 +          ff(mx,my,iz,iv,im-2)                &
                               ) * cef4(iz) * ( vp(iz,im) + 1._DP / vp(iz,im) ) &
                             + nu_s * 3._DP * ff(mx,my,iz,iv,im)                &
                             - nu_s * ksq(mx,my,iz) * Anum(ranks) * tau(ranks)  &
                               / ( Znum(ranks) * omg(iz) )**2 * ff(mx,my,iz,iv,im) & 
                               * real(iFLR, kind=DP)
                  end do
                end do
              end do
!$OMP end do nowait
            end do

          end if

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli_ct",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_LB_model


!--------------------------------------
  SUBROUTINE colli_GK_CT( ff, phi, dfdvp, im, cf )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Differential and FLR terms of test particle part in gyrokinetic collision
!                                                         with 4th order CFD
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi
    integer, intent(in) :: im
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: dfdvp

    real(kind=DP) :: cef1, cef2
    real(kind=DP), dimension(-nz:nz-1) :: cef3, cef4, cef5
    integer  ::  mx, my, iz, iv, is1


!$OMP master
                                      call clock_sta(1311)
                                    ! call fapp_start("literm_colli_ct",1311,1)
!$OMP end master


      cef1   = 1._DP / ( 12._DP * dv * dv )
      cef2   = 1._DP / ( 12._DP * dv )

!$OMP do schedule (dynamic)
      do iz = -nz, nz-1
        cef3(iz)   = 1._DP / ( 12._DP * dvp(iz) * dvp(iz) )
        cef4(iz)   = 1._DP / ( 12._DP * dvp(iz) )
        cef5(iz)   = 1._DP / ( 144._DP * dvp(iz) * dv )
      end do
!$OMP end do nowait

      if ( rankm == 0 .AND. im == 0 ) then

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  cf(mx,my,iz,iv) =                                                &
                             ( -          ff(mx,my,iz,iv-2,im)                     &
                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                 ) * 0.5_DP                                        &          
                           + ( - 30._DP * ff(mx,my,iz,iv,im  )                     &
                               + 32._DP * ff(mx,my,iz,iv,im+1)                     &
                               -  2._DP * ff(mx,my,iz,iv,im+2)                     &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                 )                                                 &           
                           + ( +          ff(mx,my,iz,iv-2,im)                     &
                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im) * 2._DP*vl(iv)**2 )           &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           - adbtc0(mx,my,iz,iv,im)*phi(mx,my,iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else if ( rankm == 0 .AND. im == 1 ) then

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  cf(mx,my,iz,iv) =                                                &
                             ( -          ff(mx,my,iz,iv-2,im)                     &
                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &          
                           + ( + 16._DP * ff(mx,my,iz,iv,im-1)                     &
                               - 31._DP * ff(mx,my,iz,iv,im  )                     &
                               + 16._DP * ff(mx,my,iz,iv,im+1)                     &
                               -          ff(mx,my,iz,iv,im+2)                     &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &           
                           + ( +          dfdvp(mx,my,iz,iv-2)                     &             
                               -  8._DP * dfdvp(mx,my,iz,iv-1)                     &
                               +  8._DP * dfdvp(mx,my,iz,iv+1)                     &
                               -          dfdvp(mx,my,iz,iv+2)                     &
                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
                               * (   nu_ps(iz,iv,im)                               & 
                                   - nu_ds(iz,iv,im) )                             &
                           + ( +          ff(mx,my,iz,iv-2,im)                     &
                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( dfdvp(mx,my,iz,iv) )                                &
                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
                                   + nu_ds(iz,iv,im)*0.5_DP                        &
                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
                                 )                                                 &
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im)                               &
                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else  ! im=[0,nm] for rankm > 0 and im=[2,nm] nm for rankm = 0  

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  cf(mx,my,iz,iv) =                                                &
                             ( -          ff(mx,my,iz,iv-2,im)                     &
                               + 16._DP * ff(mx,my,iz,iv-1,im)                     &
                               - 30._DP * ff(mx,my,iz,iv  ,im)                     &
                               + 16._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &          
                           + ( -          ff(mx,my,iz,iv,im-2)                     &
                               + 16._DP * ff(mx,my,iz,iv,im-1)                     &
                               - 30._DP * ff(mx,my,iz,iv,im  )                     &
                               + 16._DP * ff(mx,my,iz,iv,im+1)                     &
                               -          ff(mx,my,iz,iv,im+2)                     &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &           
                           + ( +          dfdvp(mx,my,iz,iv-2)                     &             
                               -  8._DP * dfdvp(mx,my,iz,iv-1)                     &
                               +  8._DP * dfdvp(mx,my,iz,iv+1)                     &
                               -          dfdvp(mx,my,iz,iv+2)                     &
                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
                               * (   nu_ps(iz,iv,im)                               & 
                                   - nu_ds(iz,iv,im) )                             &
                           + ( +          ff(mx,my,iz,iv-2,im)                     &
                               -  8._DP * ff(mx,my,iz,iv-1,im)                     &
                               +  8._DP * ff(mx,my,iz,iv+1,im)                     &
                               -          ff(mx,my,iz,iv+2,im)                     &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( dfdvp(mx,my,iz,iv) )                                &
                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
                                   + nu_ds(iz,iv,im)*0.5_DP                        &
                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
                                 )                                                 &
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im)                               &
                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
                end do
              end do
            end do
          end do
!$OMP end do nowait

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli_ct",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_GK_CT


!--------------------------------------
  SUBROUTINE colli_GK_CT6( ff, phi, dfdvp, im, cf )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Differential and FLR terms of test particle part in gyrokinetic collision
!                                                         with 6th order CFD
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)                       :: phi
    integer, intent(in) :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb) :: dfdvp

    real(kind=DP) :: cef1, cef2
    real(kind=DP), dimension(-nz:nz-1) :: cef3
    integer  ::  mx, my, iz, iv, is1


!$OMP master
                                      call clock_sta(1311)
                                    ! call fapp_start("literm_colli_ct",1311,1)
!$OMP end master


      cef1   = 1._DP / ( 90._DP * dv * dv )
      cef2   = 1._DP / ( 60._DP * dv )

      do iz = -nz, nz-1
        cef3(iz)   = 1._DP / ( 90._DP * dvp(iz) * dvp(iz) )
      end do

      if ( rankm == 0 .AND. im == 0 ) then

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                   cf(mx,my,iz,iv) =                                               &
                             ( +           ff(mx,my,iz,iv-3,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                 ) * 0.5_DP                                        &          
                           + ( - 245._DP * ff(mx,my,iz,iv,im  )                    &
                               + 270._DP * ff(mx,my,iz,iv,im+1)                    &
                               -   27_DP * ff(mx,my,iz,iv,im+2)                    &
                               +   2._DP * ff(mx,my,iz,iv,im+3)                    &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                 )                                                 &           
                           + ( -           ff(mx,my,iz,iv-3,im)                    &
                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im) * 2._DP*vl(iv)**2 )           &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           - adbtc0(mx,my,iz,iv,im)*phi(mx,my,iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else if ( rankm == 0 .AND. im == 1 ) then

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                   cf(mx,my,iz,iv) =                                               &
                             ( +           ff(mx,my,iz,iv-3,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &          
                           + ( + 135._DP  * ff(mx,my,iz,iv,im-1)                   &
                               - 258.5_DP * ff(mx,my,iz,iv,im  )                   &
                               + 136._DP  * ff(mx,my,iz,iv,im+1)                   &
                               - 13.5_DP  * ff(mx,my,iz,iv,im+2)                   &
                               +            ff(mx,my,iz,iv,im+3)                   &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &           
                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
                               +           dfdvp(mx,my,iz,iv+3)                    &
                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
                               * (   nu_ps(iz,iv,im)                               & 
                                   - nu_ds(iz,iv,im) )                             &
                           + ( -           ff(mx,my,iz,iv-3,im)                    &
                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( dfdvp(mx,my,iz,iv) )                                &
                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
                                   + nu_ds(iz,iv,im)*0.5_DP                        &
                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
                                 )                                                 &
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im)                               &
                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else if ( rankm == 0 .AND. im == 2 ) then

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                   cf(mx,my,iz,iv) =                                               &
                             ( +           ff(mx,my,iz,iv-3,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &          
                           + ( - 13.5_DP * ff(mx,my,iz,iv,im-2)                    &
                               + 136._DP * ff(mx,my,iz,iv,im-1)                    &
                               - 245._DP * ff(mx,my,iz,iv,im  )                    &
                               + 135._DP * ff(mx,my,iz,iv,im+1)                    &
                               - 13.5_DP * ff(mx,my,iz,iv,im+2)                    &
                               +           ff(mx,my,iz,iv,im+3)                    &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &           
                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
                               +           dfdvp(mx,my,iz,iv+3)                    &
                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
                               * (   nu_ps(iz,iv,im)                               & 
                                   - nu_ds(iz,iv,im) )                             &
                           + ( -           ff(mx,my,iz,iv-3,im)                    &
                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( dfdvp(mx,my,iz,iv) )                                &
                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
                                   + nu_ds(iz,iv,im)*0.5_DP                        &
                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
                                 )                                                 &
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im)                               &
                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
                end do
              end do
            end do
          end do
!$OMP end do nowait

      else  ! im=[0,nm] for rankm > 0 and im=[3,nm] for rankm = 0  

        is1 = ranks
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                   cf(mx,my,iz,iv) =                                               &
                             ( +           ff(mx,my,iz,iv-3,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv-2,im)                    &
                               + 135._DP * ff(mx,my,iz,iv-1,im)                    &
                               - 245._DP * ff(mx,my,iz,iv  ,im)                    &
                               + 135._DP * ff(mx,my,iz,iv+1,im)                    &
                               - 13.5_DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef1                                              &
                               * (   nu_ps(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ds(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &          
                           + ( +           ff(mx,my,iz,iv,im-3)                    &
                               - 13.5_DP * ff(mx,my,iz,iv,im-2)                    &
                               + 135._DP * ff(mx,my,iz,iv,im-1)                    &
                               - 245._DP * ff(mx,my,iz,iv,im  )                    &
                               + 135._DP * ff(mx,my,iz,iv,im+1)                    &
                               - 13.5_DP * ff(mx,my,iz,iv,im+2)                    &
                               +           ff(mx,my,iz,iv,im+3)                    &
                             ) * cef3(iz)                                          &
                               * (   nu_ds(iz,iv,im)*vl(iv)**2                     & 
                                   + nu_ps(iz,iv,im)*vp(iz,im)**2                  &
                                 ) * 0.5_DP                                        &           
                           + ( -           dfdvp(mx,my,iz,iv-3)                    &             
                               +   9._DP * dfdvp(mx,my,iz,iv-2)                    &
                               -  45._DP * dfdvp(mx,my,iz,iv-1)                    &
                               +  45._DP * dfdvp(mx,my,iz,iv+1)                    &
                               -   9._DP * dfdvp(mx,my,iz,iv+2)                    &
                               +           dfdvp(mx,my,iz,iv+3)                    &
                             ) * cef2 * vl(iv) * vp(iz,im)                         & 
                               * (   nu_ps(iz,iv,im)                               & 
                                   - nu_ds(iz,iv,im) )                             &
                           + ( -           ff(mx,my,iz,iv-3,im)                    &
                               +   9._DP * ff(mx,my,iz,iv-2,im)                    &
                               -  45._DP * ff(mx,my,iz,iv-1,im)                    &
                               +  45._DP * ff(mx,my,iz,iv+1,im)                    &
                               -   9._DP * ff(mx,my,iz,iv+2,im)                    &
                               +           ff(mx,my,iz,iv+3,im)                    &
                             ) * cef2 * nu_gs(iz,iv,im)*vl(iv)                     & 
                           + ( dfdvp(mx,my,iz,iv) )                                &
                               * (   nu_gs(iz,iv,im)*vp(iz,im)                     & 
                                   + nu_ds(iz,iv,im)*0.5_DP                        &
                                     * (vl(iv)**2/vp(iz,im)+vp(iz,im))             &
                                 )                                                 &
                           + ( nu_hs(iz,iv,im)*xxa(iz,iv,im)**2*2._DP              & 
                                - 0.25_DP*ksq(mx,my,iz)*Anum(is1)*tau(is1)         &
                                 * ( nu_ds(iz,iv,im)                               &
                                      * (2._DP*vl(iv)**2 + vp(iz,im)**2)           &
                                    + nu_ps(iz,iv,im)*vp(iz,im)**2 )               &
                                 / (Znum(is1)*omg(iz))**2 * real(iFLR, kind=DP)    &
                             ) * ff(mx,my,iz,iv,im)                                &
                           + adbtc(mx,my,iz,iv,im)*phi(mx,my,iz)                    
                end do
              end do
            end do
          end do
!$OMP end do nowait

      end if

!$OMP master
                                    ! call fapp_stop("literm_colli_ct",1311,1)
                                      call clock_end(1311)
!$OMP end master


  END SUBROUTINE colli_GK_CT6


!--------------------------------------
  SUBROUTINE colli_GK_DT( moment_ab_wk, im, cf )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Non-isothermal terms of test particle part and field particle part 
!                                                      in gyrokinetic collision
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ab_wk
    integer, intent(in)                                :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    integer  ::  mx, my, iz, iv, is1, is2


!$OMP master
                                      call clock_sta(1312)
                                    ! call fapp_start("literm_colli_dt",1312,1)
!$OMP end master


     if ( iFLR == 1 ) then  ! full-GK

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + x_tst(1,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(1,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(2,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(2,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + x_tst(3,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(3,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(4,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(4,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(5,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(5,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + x_tst(6,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(6,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     else if ( iFLR == 0 ) then ! DK-limit

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + x_tst(1,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(1,mx,my,iz,is2)      &
                                    + x_tst(3,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(3,mx,my,iz,is2)      &
                                    + x_tst(4,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(4,mx,my,iz,is2)      &
                                    + x_tst(6,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(6,mx,my,iz,is2)          
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     end if


!$OMP master
                                    ! call fapp_stop("literm_colli_dt",1312,1)
                                      call clock_end(1312)
!$OMP end master


  END SUBROUTINE colli_GK_DT


!--------------------------------------
  SUBROUTINE colli_GK_CF_DT(moment_ba_wk, moment_ab_wk, im, cf)
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Field particle and non-isothermal parts in gyrokinetic collision
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ba_wk, moment_ab_wk
    integer, intent(in)                                :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    integer  ::  mx, my, iz, iv, is1, is2


!$OMP master
                                      call clock_sta(1313)
                                    ! call fapp_start("literm_colli_cf",1313,1)
!$OMP end master

     if ( iFLR == 1 ) then  ! full-GK

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + y_fld(1,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(1,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(2,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(2,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + y_fld(3,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(3,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(4,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(4,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(5,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(5,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + y_fld(6,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(6,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  &

                                    + x_tst(1,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(1,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(2,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(2,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + x_tst(3,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(3,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(4,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(4,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + x_tst(5,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(5,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + x_tst(6,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(6,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     else if ( iFLR == 0 ) then ! DK-limit

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + y_fld(1,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(1,mx,my,iz,is2)      &
                                    + y_fld(3,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(3,mx,my,iz,is2)      &
                                    + y_fld(4,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(4,mx,my,iz,is2)      &
                                    + y_fld(6,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(6,mx,my,iz,is2)      &

                                    + x_tst(1,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(1,mx,my,iz,is2)      &
                                    + x_tst(3,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(3,mx,my,iz,is2)      &
                                    + x_tst(4,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(4,mx,my,iz,is2)      &
                                    + x_tst(6,iz,iv,im,is1,is2)          &
                                     * moment_ab_wk(6,mx,my,iz,is2)          
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     end if

!$OMP master
                                    ! call fapp_stop("literm_colli_cf",1313,1)
                                      call clock_end(1313)
!$OMP end master


  END SUBROUTINE colli_GK_CF_DT


!--------------------------------------
  SUBROUTINE colli_GK_CF(moment_ba_wk, im, cf)
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Field particle part in gyrokinetic collision
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(1:6,-nx:nx,0:ny,-nz:nz-1,0:ns-1)       :: moment_ba_wk
    integer, intent(in)                                :: im
    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv)           :: cf

    integer  ::  mx, my, iz, iv, is1, is2


!$OMP master
                                      call clock_sta(1313)
                                    ! call fapp_start("literm_colli_cf",1313,1)
!$OMP end master

     if ( iFLR == 1 ) then  ! full-GK

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + y_fld(1,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(1,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(2,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(2,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + y_fld(3,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(3,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(4,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(4,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  & 
                                    + y_fld(5,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(5,mx,my,iz,is2)      &
                                      * j1(mx,my,iz,im)                  & 
                                    + y_fld(6,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(6,mx,my,iz,is2)      &
                                      * j0(mx,my,iz,im)                  
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     else if ( iFLR == 0 ) then ! DK-limit

       is1 = ranks
         do is2 = 0, ns-1
!%%%% do schedule (dynamic)
!$OMP do
           do iv = 1, 2*nv
             do iz = -nz, nz-1
               do my = ist_y, iend_y
                 do mx = -nx, nx
                   cf(mx,my,iz,iv) = cf(mx,my,iz,iv)                     &
                                    + y_fld(1,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(1,mx,my,iz,is2)      &
                                    + y_fld(3,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(3,mx,my,iz,is2)      &
                                    + y_fld(4,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(4,mx,my,iz,is2)      &
                                    + y_fld(6,iz,iv,im,is1,is2)          &
                                     * moment_ba_wk(6,mx,my,iz,is2)         
                end do
              end do
            end do
          end do
!%%%% end do nowait
!$OMP end do
        end do

     end if

!$OMP master
                                    ! call fapp_stop("literm_colli_cf",1313,1)
                                      call clock_end(1313)
!$OMP end master


  END SUBROUTINE colli_GK_CF


!--------------------------------------
  SUBROUTINE colli_moment_calc( hh, phi, ww )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Moment calculations for gyrokinetic collision: local velocity moment part 
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: phi

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww 

    real(kind=DP) :: v2a, v2b, dflg
    complex(kind=DP) :: wf1, wf2
    integer :: mx, my, iz, iv, im, is1, is2, ii

!$OMP master
                                      call clock_sta(1314)
                                    ! call fapp_start("literm_colli_mom",1314,1)
!$OMP end master

      if ( rankm == 0 ) then

!$OMP do collapse(2) schedule(dynamic)
      do ii = 1, 6
        do is2 = 0, ns-1
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*vfunc(mx,my,iz,iv,im,is2,ii)
                  end do
                end do
              end do
            end do
          end do

! for edge compensation
!         im  = 1
          do im = 1, 1
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    v2a = vfunc(mx,my,iz,iv,im,is2,ii)
                    v2b = vfunc(mx,my,iz,iv,im+1,is2,ii)
                    wf1 = hh(mx,my,iz,iv,im)   
                    wf2 = hh(mx,my,iz,iv,im+1)
                    ww(mx,my,iz,is2,ii)  = ww(mx,my,iz,is2,ii)            &
                          - ( - wf1/12._DP*v2a + ( wf2*v2b - wf1*2._DP*v2a )*11._DP/720._DP ) 
                  end do
                end do
              end do
            end do
          end do
        end do
      enddo
!$OMP enddo nowait

      else

!$OMP do collapse(2) schedule(dynamic)
      do ii=1,6
        do is2 = 0, ns-1
          do im = 0, nm
            do iv = 1, 2*nv
              do iz = -nz, nz-1
                do my = ist_y, iend_y
                  do mx = -nx, nx
                    ww(mx,my,iz,is2,ii) = ww(mx,my,iz,is2,ii) + hh(mx,my,iz,iv,im)*vfunc(mx,my,iz,iv,im,is2,ii)
                  end do
                end do
              end do
            end do
          end do
        end do
      enddo
!$OMP enddo nowait

      end if


!$OMP master
                                    ! call fapp_stop("literm_colli_mom",1314,1)
                                      call clock_end(1314)
!$OMP end master

  END SUBROUTINE colli_moment_calc


!--------------------------------------
  SUBROUTINE colli_moment_redc( ww, wn )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Moment calculations for gyrokinetic collision: All_reduce_part
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww
 
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: wn


!$OMP master
                                           call clock_sta(1315)
                                         ! call fapp_start("literm_colli_ar",1315,1)
!$OMP end master

      call MPI_Allreduce( ww, wn, nxyz*ns*6, MPI_DOUBLE_COMPLEX, &
                          MPI_SUM, vel_comm_world, ierr_mpi )

!$OMP master
                                         ! call fapp_stop("literm_colli_ar",1315,1)
                                           call clock_end(1315)
!$OMP end master

  END SUBROUTINE colli_moment_redc


!--------------------------------------
  SUBROUTINE colli_comm_alltoall( wm, wn )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    Inter-species communication of moment quantities for field particle part
!       with MPI_AlltoAll
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)  :: wm
    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)  :: wn

    complex(kind=DP),              & 
      dimension(-nx:nx,0:ny,-nz:nz-1,1:6,0:ns-1)  :: send_buff, recv_buff

    integer :: mx, my, iz, is, ii
    integer :: datasize, datasize_ns


!$OMP master
                                       call clock_sta(1316)
                                     ! call fapp_start("literm_colli_com",1316,1)
!$OMP end master


      datasize = (2*nx+1)*(ny+1)*(2*nz)*6
      datasize_ns = (2*nx+1)*(ny+1)*(2*nz)*6*ns

      if ( vel_rank == 0 ) then

        do ii = 1, 6
          do is = 0, ns-1
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
!                  send_buff(mx,my,iz,ii,is)  = real(ranks,kind=DP) + ii ! for debug
                  send_buff(mx,my,iz,ii,is)  = wm(mx,my,iz,is,ii)
                end do
              end do 
            end do
          end do
        end do


          call MPI_Alltoall( send_buff(-nx,ist_y,-nz,1,0), datasize, MPI_DOUBLE_COMPLEX,  &
                             recv_buff(-nx,ist_y,-nz,1,0), datasize, MPI_DOUBLE_COMPLEX,  & 
                             col_comm_world, ierr_mpi  )
 
!! --- for debug 
!        write(unit=8000+ranks,fmt="(I3,SP,6ES24.15e3)") ranks, real(recv_buff(0,1,1,1,0)), real(recv_buff(0,1,1,6,0)), &
!                                                           real(recv_buff(0,1,1,1,1)), real(recv_buff(0,1,1,6,1)), &
!                                                           real(recv_buff(0,1,1,1,2)), real(recv_buff(0,1,1,6,2))
!        write(unit=8000+ranks,fmt=*)

      end if


      call MPI_Bcast( recv_buff(-nx,ist_y,-nz,1,0), datasize_ns, MPI_DOUBLE_COMPLEX, & 
                      0, vel_comm_world, ierr_mpi  ) 


      do is = 0, ns-1
        do ii = 1, 6
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                wn(mx,my,iz,is,ii) = recv_buff(mx,my,iz,ii,is)
              end do
            end do 
          end do
        end do         
      end do         

!$OMP master
                                     ! call fapp_stop("literm_colli_com",1316,1)
                                       call clock_end(1316)
!$OMP end master

  END SUBROUTINE colli_comm_alltoall


!--------------------------------------
  SUBROUTINE colli_dfdvp( ff, dfdvp )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    calculation of df/dv_perp term with 4th order CFD
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: dfdvp

    real(kind=DP), dimension(-nz:nz-1) :: cef4
    integer  ::  mx, my, iz, iv, im


!$OMP master
                                      call clock_sta(1317)
                                    ! call fapp_start("literm_colli_dvp",1317,1)
!$OMP end master


      do iz = -nz, nz-1
        cef4(iz)   = 1._DP / ( 12._DP * dvp(iz) )
      end do


      if ( rankm == 0 ) then 

        im = 0
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) = (0._DP, 0._DP) 
                end do
              end do
            end do
          end do
!$OMP end do nowait

        im = 1
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( -  8._DP * ff(mx,my,iz,iv,im-1)                &
                               +          ff(mx,my,iz,iv,im  )                &
                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
                               -          ff(mx,my,iz,iv,im+2)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait


        do im = 2, nm
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( +          ff(mx,my,iz,iv,im-2)                &
                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
                               -          ff(mx,my,iz,iv,im+2)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do

      else   

        do im = 0, nm
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( +          ff(mx,my,iz,iv,im-2)                &
                               -  8._DP * ff(mx,my,iz,iv,im-1)                &
                               +  8._DP * ff(mx,my,iz,iv,im+1)                &
                               -          ff(mx,my,iz,iv,im+2)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait
        end do
     
      end if 


!$OMP master
                                    ! call fapp_stop("literm_colli_dvp",1317,1)
                                      call clock_end(1317)
!$OMP end master


  END SUBROUTINE colli_dfdvp


!--------------------------------------
  SUBROUTINE colli_dfdvp6( ff, dfdvp )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    calculation of df/dv_perp term with 6th order CFD
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

    complex(kind=DP), intent(out), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: dfdvp

    real(kind=DP), dimension(-nz:nz-1) :: cef4
    integer  ::  mx, my, iz, iv, im


!$OMP master
                                      call clock_sta(1317)
                                    ! call fapp_start("literm_colli_dvp",1317,1)
!$OMP end master


      do iz = -nz, nz-1
        cef4(iz)   = 1._DP / ( 60._DP * dvp(iz) )
      end do


      if ( rankm == 0 ) then 

        im = 0
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) = (0._DP, 0._DP) 
                end do
              end do
            end do
          end do
!$OMP end do nowait

        im = 1
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( - 45._DP * ff(mx,my,iz,iv,im-1)                &
                               +  9._DP * ff(mx,my,iz,iv,im  )                &
                               + 44._DP * ff(mx,my,iz,iv,im+1)                &
                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
                               +          ff(mx,my,iz,iv,im+3)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait

        im = 2
!$OMP do schedule (dynamic)
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( +  9._DP * ff(mx,my,iz,iv,im-2)                &
                               - 46._DP * ff(mx,my,iz,iv,im-1)                &
                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
                               +          ff(mx,my,iz,iv,im+3)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
!$OMP end do nowait


!$OMP do schedule (dynamic)
        do im = 3, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( -          ff(mx,my,iz,iv,im-3)                &
                               +  9._DP * ff(mx,my,iz,iv,im-2)                &
                               - 45._DP * ff(mx,my,iz,iv,im-1)                &
                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
                               +          ff(mx,my,iz,iv,im+3)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
        end do
!$OMP end do nowait

      else   

!$OMP do schedule (dynamic)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  dfdvp(mx,my,iz,iv,im) =                                     &
                             ( -          ff(mx,my,iz,iv,im-3)                &
                               +  9._DP * ff(mx,my,iz,iv,im-2)                &
                               - 45._DP * ff(mx,my,iz,iv,im-1)                &
                               + 45._DP * ff(mx,my,iz,iv,im+1)                &
                               -  9._DP * ff(mx,my,iz,iv,im+2)                &
                               +          ff(mx,my,iz,iv,im+3)                &
                             ) * cef4(iz)
                end do
              end do
            end do
          end do
        end do
!$OMP end do nowait
     
      end if 


!$OMP master
                                    ! call fapp_stop("literm_colli_dvp",1317,1)
                                      call clock_end(1317)
!$OMP end master


  END SUBROUTINE colli_dfdvp6


!--------------------------------------
  SUBROUTINE colli_zeroset( cff )
!--------------------------------------
!-------------------------------------------------------------------------------
!
!    zero clear for collision terms 
!
!    by M. Nakata, M. Nunami, April 2014
!
!-------------------------------------------------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm)    :: cff

    integer  ::  mx, my, iz, iv, im

!$OMP master
                                      call clock_sta(1318)
                                    ! call fapp_start("literm_colli_0",1318,1)
!$OMP end master

!$OMP do schedule (dynamic)
        do im = 0, nm
          do iv = 1, 2*nv
            do iz = -nz, nz-1
              do my = ist_y, iend_y
                do mx = -nx, nx
                  cff(mx,my,iz,iv,im) = ( 0._DP, 0._DP )     
                end do
              end do
            end do
          end do
        end do
!$OMP end do nowait 


!$OMP master
                                    ! call fapp_stop("literm_colli_0",1318,1)
                                      call clock_end(1318)
!$OMP end master


  END SUBROUTINE colli_zeroset


!--------------------------------------
  SUBROUTINE colli_hhset(hh,phi,ff)
!--------------------------------------
    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz:nz-1)               :: phi
    complex(kind=DP), intent(in), &
      dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff
    real(kind=DP) :: dflg
    integer :: mx, my, iz, iv, im, is1

!$OMP master
                                      call clock_sta(1319)
                                    ! call fapp_start("literm_colli_hwset",1319,1)
!$OMP end master

      dflg = real(1-icheck,kind=DP) 

      is1 = ranks
!$OMP do collapse(2)
      do im = 0, nm
        do iv = 1, 2*nv
          do iz = -nz, nz-1
            do my = ist_y, iend_y
              do mx = -nx, nx
                hh(mx,my,iz,iv,im) = ( ff(mx,my,iz,iv,im) + dflg*sgn(is1)*j0(mx,my,iz,im)*phi(mx,my,iz)   &
                                                                * fmx(iz,iv,im)*Znum(is1)/tau(is1) )      &
                                    * vp(iz,im) * dvp(iz) * dv * twopi
              end do
            end do
          end do
        end do
      end do
!$OMP enddo nowait

!$OMP master
                                    ! call fapp_stop("literm_colli_hwset",1319,1)
                                      call clock_end(1319)
!$OMP end master


  END SUBROUTINE colli_hhset


!--------------------------------------
  SUBROUTINE colli_wwset(ww)
!--------------------------------------
    complex(kind=DP), &
      dimension(-nx:nx,0:ny,-nz:nz-1,0:ns-1,1:6)        :: ww 
    integer :: mx, my, iz, is1, ii

!$OMP master
                                      call clock_sta(1319)
                                    ! call fapp_start("literm_colli_hwset",1319,1)
!$OMP end master

!$OMP do collapse(2)
    do ii = 1, 6
      do is1 = 0, ns-1
        do iz = -nz, nz-1
          do my = ist_y, iend_y
            do mx = -nx, nx
              ww(mx,my,iz,is1,ii) = ( 0._DP, 0._DP )
            end do 
          end do 
        end do 
      end do 
    end do 
!$OMP enddo

!$OMP master
                                    ! call fapp_stop("literm_colli_hwset",1319,1)
                                      call clock_end(1319)
!$OMP end master


  END SUBROUTINE colli_wwset


END MODULE GKV_colli
