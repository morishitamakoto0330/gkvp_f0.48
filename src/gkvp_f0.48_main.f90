PROGRAM GKV_main
!-------------------------------------------------------------------------------
!
!    Nonlinear gyrokinetic Vlasov code in a flux tube geometry
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!               f0.30 ( S. Maeyama, March 2013 )
!                     updated for electromagnetic, multi-species,
!                                 MHD equilibrium, 5D-parallelization
!               f0.40 ( M. Nakata, June 2014 )
!                     updated for realistic tokamak equilibrium, 
!                                 multi-species collision 
!
!      Hierarchy of the modules (The lower should be complied earlier)
!
!        main
!         |
!        set
!         |
!        advance, dtc, out
!         |
!        exb, trans
!         |
!        bndry, colli, fft, fld, freq, zfilter
!         |
!        clock, intgrl, tips, vmecin, igs
!         |
!        header, mpienv, math
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv
  use GKV_set,   only: set_init, set_close
  use GKV_clock, only: clock_timer, clock_sta, clock_end, clock_reset
  use GKV_out,   only: out_cntrl, contnu
  use GKV_dtc,   only: dtc_cntrl
  use GKV_fld,   only: fld_esfield
  use GKV_advnc, only: advnc_rkgsteps
  use GKV_fft,   only: fft_pre
  use GKV_freq,  only: freq_conv
  use mpi !fj

  implicit none

  complex(kind=DP), &
    dimension(-nx:nx,0:ny,-nz-nzb:nz-1+nzb,1-nvb:2*nv+nvb,0-nvb:nm+nvb) :: ff

  complex(kind=DP), dimension(-nx:nx,0:ny,-nz:nz-1)       :: Al, phi

  complex(kind=DP), &
    dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: hh, dh, cft, cff, ef

  real(kind=DP) :: time

  integer :: loop, iflg

  integer :: ny_size, nxw_sz, nwk


    call mpienv_init( nprocw, nprocz, nprocv, nprocm, nprocs )

! ---- set y range --------------------------
    ny_size = global_ny + 1 
    if( mod(ny_size,nprocw) == 0 )  then
      nwk    = ny_size / nprocw
    else
      nwk    = ny_size / nprocw + 1
    endif
    !--- global index range ---------------- 
    ist_y_g  = nwk*rankw
    iend_y_g = min( nwk*(rankw+1)-1, (ny_size-1) )
    nsize_y  = iend_y_g - ist_y_g + 1
    !--- local index range ---------------- 
    ist_y    = 0
    iend_y   = iend_y_g - ist_y_g

    if( rankw == 0 )   then
       ist1_y    = 1
    else 
       ist1_y    = 0
    endif

! ---- set xw range ---------------------
    nxw_sz = 2*nxw
    if( mod(nxw_sz,nprocw) == 0 )  then
      nwk    = nxw_sz / nprocw
    else
      nwk    = nxw_sz / nprocw + 1
    endif
    !--- global index range ----------------
    ist_xw_g  = nwk*rankw
    iend_xw_g = min( nwk*(rankw+1)-1, (nxw_sz-1) )
    nsize_xw  = iend_xw_g - ist_xw_g + 1
    !--- local index range ----------------
    ist_xw    = 0
    iend_xw   = iend_xw_g - ist_xw_g

    call clock_timer( 0, iflg )
                                           call clock_sta(1)
                                         ! call fapp_start("pre",1,1)
    call fft_pre( )
    call set_init( ff, phi, Al, hh, dh, cft, cff, ef, time )
      write( olog, * ) " # simulation is started at t = ", time

    call out_cntrl( ff, phi, Al, dh, cft, cff, time, 0 )

    if ( adapt_dt ) call dtc_cntrl( time, 0 )
                                         ! call fapp_stop("pre",1,1)
                                           call clock_end(1)
                                           call clock_reset
    
    loop   = 0
                                           call clock_sta(2)
                                         ! call fipp_start
                                         ! call fapp_start("timesteploop",2,1)
    do

      if ( time > tend - eps ) exit

      time   = time + dt
      loop   = loop + 1

      call advnc_rkgsteps( ff, phi, Al, hh, dh, cft, cff, ef )
      
                                           call clock_sta(10)
                                         ! call fapp_start("output",10,1)
      call out_cntrl( ff, phi, Al, dh, cft, cff, time, 1 )
      if ( adapt_dt ) call dtc_cntrl( time, 1 )
      if ( calc_type == "lin_freq" .and. all(freq_conv) ) then
        write( olog, * ) " # Growth rate and frequency are well converged."
        exit
      end if
                                         ! call fapp_stop("output",10,1)
                                           call clock_end(10)

! --- output continu file every 10000 steps
      if (mod(loop+10000,10000) == 0 ) then 
                                           call clock_sta(16)
        write( olog, * ) "# check-point at time = ", time
        call flush(olog)
        call contnu ( ff, time )
        call flush(ocnt)
        call flush(ofxv)
        if ( vel_rank == 0 ) then
          call flush(omom)
        end if
        if ( ranks == 0 .AND. vel_rank == 0 ) then
          call flush(ophi)
          call flush(oAl)
        end if
        if ( zsp_rank == 0 .AND. vel_rank == 0 ) then
          call flush(otrn)
        end if
        if( rankg == 0 ) then
          call flush(odtc)
          call flush(oeng)
          call flush(omen)
          call flush(owes)
          call flush(owem)
          if ( trim(calc_type) == "linear" .or. &
               trim(calc_type) == "lin_freq" ) then
            call flush(ofrq)
          end if
        end if
        if( rank == 0 ) then
          call flush(obln)
          call flush(oges)
          call flush(ogem)
          call flush(oqes)
          call flush(oqem)
        end if
                                           call clock_end(16)
      end if
! ---
      call clock_timer( 1, iflg )
      
      if( iflg == 1 ) exit

    end do
                                         ! call fapp_stop("timesteploop",2,1)
                                         ! call fipp_stop
                                           call clock_end(2)

                                           call clock_sta(3)
                                         ! call fapp_start("post",3,1)
    call out_cntrl( ff, phi, Al, dh, cft, cff, time, 2 )
      write( olog, * ) " # simulation is stopped at t = ", time
                                         ! call fapp_stop("post",3,1)
                                           call clock_end(3)
    call clock_timer( 2, iflg )

    call set_close

    call MPI_Finalize ( ierr_mpi )

  stop


END PROGRAM GKV_main
