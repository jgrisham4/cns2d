!===============================================================================
! This module contains different limiters for use in the
! FVM solver.  Because the reconstruction is done in an
! unstructured manner, only the barth limiter is used.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

module limiters
  use solver_class,    only : solver
  use ieee_arithmetic, only : ieee_is_finite
  implicit none
  private
  public :: minmod, vanleer, barth

  contains

    !---------------------------------------------------------------------------
    ! minmod limiter
    !---------------------------------------------------------------------------
    double precision function minmod(r) result(phi)
      implicit none
      double precision, intent(in) :: r
      if (r.le.0.0d0) then
        phi = 0.0d0
      else if (r.lt.1.0d0) then
        phi = r
      else
        phi = 1.0d0
      end if
    end function minmod

    !---------------------------------------------------------------------------
    ! Another version of the minmod limiter
    !---------------------------------------------------------------------------
    double precision function minmod2(a,b) result(sigma)
      implicit none
      double precision, intent(in) :: a,b
      if ((abs(a).gt.abs(b)).and.(a*b.gt.0.0d0)) then
        sigma = a
      else if ((abs(b).gt.abs(a)).and.(a*b.gt.0.0d0)) then
        sigma = b
      else
        sigma = 0.0d0
      end if
    end function minmod2

    !---------------------------------------------------------------------------
    ! van Leer limiter
    !---------------------------------------------------------------------------
    double precision function vanleer(r) result(phi)
      implicit none
      double precision, intent(in) :: r
      if (r.le.0.0d0) then
        phi = 0.0d0
      else if (isnan(r).or..not.ieee_is_finite(r)) then
        phi = 2.0d0
      else
        phi = 2.0d0*r/(r+1.0d0)
      end if
    end function vanleer

    !---------------------------------------------------------------------------
    ! Subroutine for the Barth-Jespersen limiter
    !
    ! NOTES:
    ! - This function assumes that the gradient of the
    !   conserved variables has already been computed.
    !---------------------------------------------------------------------------
    pure function barth(s,i,j,gradu,umax,umin) result(ph)
      implicit none
      type(solver), intent(in)                  :: s
      integer, intent(in)                       :: i,j
      double precision, allocatable, intent(in) :: gradu(:,:,:),umax(:,:,:),umin(:,:,:)
      double precision                          :: ph(4)
      double precision                          :: rvec(2),unodes(4,4),phibar(4,4)
      integer                                   :: k,l

      ! Node numbering:
      !
      !  4           3
      !   o---------o
      !   |         |
      !   |    +    |
      !   |         |
      !   o---------o
      !  1           2

      ! Extrapolating state to vertices of element
      ! Barth's paper says that the extrema occur at the vertices.
      ! I agree, but I'm not reconstructing the state at the vertices.
      ! Blazek's book says to use the state extrapolated to the face midpoints
      ! for the cell-centered scheme.
      ! Point 1
      !rvec(1) = s%grid%x(i,j) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i,j) - s%grid%elem(i,j)%yc
      ! Left side
      rvec(1) = s%grid%edges_h(i,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_h(i,j)%ym - s%grid%elem(i,j)%xc
      unodes(:,1) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 2
      !rvec(1) = s%grid%x(i+1,j) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i+1,j) - s%grid%elem(i,j)%yc
      ! Right side
      rvec(1) = s%grid%edges_h(i+1,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_h(i+1,j)%ym - s%grid%elem(i,j)%yc
      unodes(:,2) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 3
      !rvec(1) = s%grid%x(i+1,j+1) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i+1,j+1) - s%grid%elem(i,j)%yc
      ! Bottom
      rvec(1) = s%grid%edges_v(i,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_v(i,j)%ym - s%grid%elem(i,j)%yc
      unodes(:,3) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 4
      !rvec(1) = s%grid%x(i,j+1) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i,j+1) - s%grid%elem(i,j)%yc
      ! Top
      rvec(1) = s%grid%edges_v(i,j+1)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_v(i,j+1)%ym - s%grid%elem(i,j)%yc
      unodes(:,4) = s%grid%elem(i,j)%u + gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Finding value of limiter at each vertex -- not anymore
      ! Finding value of limiter at each face midpoint
      ! k - vertex, l - variables
      do k=1,4
        do l=1,4
          if (unodes(l,k)-s%grid%elem(i,j)%u(l).gt.0.0d0) then
            phibar(l,k) = min(1.0d0,(umax(i,j,l)-s%grid%elem(i,j)%u(l))/(unodes(l,k)-s%grid%elem(i,j)%u(l)))
          else if (unodes(l,k)-s%grid%elem(i,j)%u(l).lt.0.0d0) then
            phibar(l,k) = min(1.0d0,(umin(i,j,l)-s%grid%elem(i,j)%u(l))/(unodes(l,k)-s%grid%elem(i,j)%u(l)))
          else
            phibar(l,k) = 1.0d0
          end if
        end do
      end do

      ! Finding the final value of the limiter for each conserved variable
      do k=1,4
        ph(k) = minval(phibar(k,:))
      end do

    end function barth

    !---------------------------------------------------------------------------
    ! Subroutine for the Venkatakrishnan's limiter
    !---------------------------------------------------------------------------
    pure function venkatakrishnan(s,i,j,gradu,umax,umin,kval) result(ph)
      implicit none
      type(solver), intent(in)                  :: s
      integer, intent(in)                       :: i,j
      double precision, allocatable, intent(in) :: gradu(:,:,:),umax(:,:,:),umin(:,:,:)
      double precision, intent(in)              :: kval
      double precision                          :: ph(4),d1max,d1min
      double precision                          :: rvec(2),d2(4,4),phibar(4,4)
      double precision                          :: eps2
      integer                                   :: k,l

      ! Node numbering:
      !
      !  4           3
      !   o---------o
      !   |         |
      !   |    +    |
      !   |         |
      !   o---------o
      !  1           2

      ! Computing parameter epsilon^2
      eps2 = (kval*sqrt(s%grid%elem(i,j)%area))**3

      ! Extrapolating state to vertices of element
      ! Barth's paper says that the extrema occur at the vertices.
      ! I agree, but I'm not reconstructing the state at the vertices.
      ! Blazek's book says to use the state extrapolated to the face midpoints
      ! for the cell-centered scheme.
      ! Point 1
      !rvec(1) = s%grid%x(i,j) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i,j) - s%grid%elem(i,j)%yc
      ! Left side
      rvec(1) = s%grid%edges_h(i,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_h(i,j)%ym - s%grid%elem(i,j)%xc
      d2(:,1) = gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 2
      !rvec(1) = s%grid%x(i+1,j) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i+1,j) - s%grid%elem(i,j)%yc
      ! Right side
      rvec(1) = s%grid%edges_h(i+1,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_h(i+1,j)%ym - s%grid%elem(i,j)%yc
      d2(:,2) = gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 3
      !rvec(1) = s%grid%x(i+1,j+1) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i+1,j+1) - s%grid%elem(i,j)%yc
      ! Bottom
      rvec(1) = s%grid%edges_v(i,j)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_v(i,j)%ym - s%grid%elem(i,j)%yc
      d2(:,3) = gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Point 4
      !rvec(1) = s%grid%x(i,j+1) - s%grid%elem(i,j)%xc
      !rvec(2) = s%grid%y(i,j+1) - s%grid%elem(i,j)%yc
      ! Top
      rvec(1) = s%grid%edges_v(i,j+1)%xm - s%grid%elem(i,j)%xc
      rvec(2) = s%grid%edges_v(i,j+1)%ym - s%grid%elem(i,j)%yc
      d2(:,4) = gradU(i,j,1:4)*rvec(1) + gradU(i,j,5:8)*rvec(2)

      ! Finding value of limiter at each vertex -- not anymore
      ! k - vertex, l - variables
      ! Finding value of limiter at each face midpoint
      ! k - face, l - variables
      do k=1,4
        do l=1,4

          ! Computing Delta_{1,max} and Delta_{1,min}
          d1min = umin(i,j,l) - s%grid%elem(i,j)%u(l)
          d1max = umax(i,j,l) - s%grid%elem(i,j)%u(l)

          ! Logic to determine limiter
          if (d2(l,k)>0.0d0) then
            phibar(l,k) = 1.0d0/d2(l,k)*(((d1max**2+eps2)*d2(l,k)+2.0d0*d2(l,k)**2*d1max)/ &
              (d1max**2 + 2.0d0*d2(l,k)**2 + d1max*d2(l,k) + eps2))
          else if (d2(l,k)<0.0d0) then
            phibar(l,k) = 1.0d0/d2(l,k)*(((d1min**2+eps2)*d2(l,k)+2.0d0*d2(l,k)**2*d1min)/ &
              (d1min**2 + 2.0d0*d2(l,k)**2 + d1min*d2(l,k) + eps2))
          else
            phibar(l,k) = 1.0d0
          end if
        end do
      end do

      ! Finding the final value of the limiter for each conserved variable
      do k=1,4
        ph(k) = minval(phibar(k,:))
      end do

    end function venkatakrishnan

end module limiters
