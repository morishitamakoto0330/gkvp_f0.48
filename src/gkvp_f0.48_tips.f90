MODULE GKV_tips
!-------------------------------------------------------------------------------
!
!    Some useful tools and tips
!
!      GKV-plus r0.3 ( T.-H.Watanabe, Jun 2011)
!
!-------------------------------------------------------------------------------

  use GKV_header
  use GKV_mpienv

  implicit none

  private

  public   tips_reality


CONTAINS


!--------------------------------------
  SUBROUTINE tips_reality( wrk )
!--------------------------------------

    complex(kind=DP), intent(inout), &
      dimension(-nx:nx,0:ny,-nz:nz-1,1:2*nv,0:nm) :: wrk

    integer :: mx

      if( rankw == 0 )  then
        do mx = 0, nx
          wrk(-mx,0,:,:,:) = conjg( wrk(mx,0,:,:,:) )
        end do
      endif


  END SUBROUTINE tips_reality


END MODULE GKV_tips
