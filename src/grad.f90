!===============================================================================
! This module contains a subroutine which can be used to
! compute the gradient of a variable on a curvilinear
! mesh.  It does so making use of the metrics of the
! general transformation.
!
! Author: James Grisham
! Date  : 12-03-2016
!===============================================================================

module grad
  use mesh_class, only : mesh
  implicit none
  private
  public :: compute_gradient

  contains

  !-----------------------------------------------------------------------------
  ! Subroutine for computing the gradient of an input
  ! array named var which has four independent variables
  !
  ! Note: This subroutine assumes that memory has been
  !       allocated for gradVar prior to this subroutine
  !       being called.
  !-----------------------------------------------------------------------------
  subroutine compute_gradient(grid,var,gradVar)
    implicit none
    type(mesh), intent(in)                       :: grid
    double precision, allocatable, intent(in)    :: var(:,:,:)
    double precision, allocatable, intent(inout) :: gradVar(:,:,:)
    double precision, dimension(4)               :: dvardxi,dvardeta
    integer                                      :: i,j

    ! Computing the gradient on the interior cells
    do j=2,grid%nelemj-1
      do i=2,grid%nelemi-1

        ! Computing the derivative in the computational domain
        dvardxi  = 0.5d0*(var(i+1,j,:)-var(i-1,j,:))
        dvardeta = 0.5d0*(var(i,j+1,:)-var(i,j-1,:))

        ! Finding the gradient in the physical domain
        gradVar(i,j,1:4) = 1.0d0/grid%elem(i,j)%detJ*(dvardxi*grid%elem(i,j)%dydeta-dvardeta*grid%elem(i,j)%dydxi)
        gradVar(i,j,5:8) = 1.0d0/grid%elem(i,j)%detJ*(dvardeta*grid%elem(i,j)%dxdxi-dvardxi*grid%elem(i,j)%dxdeta)

      end do
    end do

    ! Computing gradients on the lower horizontal boundary
    j = 1
    do i=1,grid%nelemi

      dvardeta = 0.5d0*(-3.0d0*var(i,j,:) + 4.0d0*var(i,j+1,:) - var(i,j+2,:))

      if (i.ne.1.and.i.ne.grid%nelemi) then
        dvardxi = 0.5d0*(var(i+1,j,:) - var(i-1,j,:))
      else if (i.eq.1) then
        dvardxi = 0.5d0*(-3.0d0*var(i,j,:) + 4.0d0*var(i+1,j,:) - var(i+2,j,:))
      else
        dvardxi = 0.5d0*(3.0d0*var(i,j,:) - 4.0d0*var(i-1,j,:) + var(i-2,j,:))
      end if

      ! Finding the gradient in the physical domain (lower)
      gradVar(i,j,1:4) = 1.0d0/grid%elem(i,j)%detJ*(dvardxi*grid%elem(i,j)%dydeta-dvardeta*grid%elem(i,j)%dydxi)
      gradVar(i,j,5:8) = 1.0d0/grid%elem(i,j)%detJ*(dvardeta*grid%elem(i,j)%dxdxi-dvardxi*grid%elem(i,j)%dxdeta)

    end do

    ! Computing gradients on upper horizontal boundary
    j = grid%nelemj
    do i=1,grid%nelemi

      dvardeta = 0.5d0*(3.0d0*var(i,j,:) - 4.0d0*var(i,j-1,:) + var(i,j-2,:))

      if (i.ne.1.and.i.ne.grid%nelemi) then
        dvardxi = 0.5d0*(var(i+1,j,:) - var(i-1,j,:))
      else if (i.eq.1) then
        dvardxi = 0.5d0*(-3.0d0*var(i,j,:) + 4.0d0*var(i+1,j,:) - var(i+2,j,:))
      else
        dvardxi = 0.5d0*(3.0d0*var(i,j,:) - 4.0d0*var(i-1,j,:) + var(i-2,j,:))
      end if

      ! Finding the gradient in the physical domain (lower)
      gradVar(i,j,1:4) = 1.0d0/grid%elem(i,j)%detJ*(dvardxi*grid%elem(i,j)%dydeta-dvardeta*grid%elem(i,j)%dydxi)
      gradVar(i,j,5:8) = 1.0d0/grid%elem(i,j)%detJ*(dvardeta*grid%elem(i,j)%dxdxi-dvardxi*grid%elem(i,j)%dxdeta)

    end do

    ! Computing gradients on the vertical boundaries
    do j=2,grid%nelemj-1

      ! Left boundary
      i = 1
      dvardxi  = 0.5d0*(-3.0d0*var(i,j,:) + 4.0d0*var(i+1,j,:) - var(i+2,j,:))
      dvardeta = 0.5d0*(var(i,j+1,:) - var(i,j-1,:))
      gradVar(i,j,1:4) = 1.0d0/grid%elem(i,j)%detJ*(dvardxi*grid%elem(i,j)%dydeta-dvardeta*grid%elem(i,j)%dydxi)
      gradVar(i,j,5:8) = 1.0d0/grid%elem(i,j)%detJ*(dvardeta*grid%elem(i,j)%dxdxi-dvardxi*grid%elem(i,j)%dxdeta)

      ! Right boundary
      i = grid%nelemi
      dvardxi  = 0.5d0*(3.0d0*var(i,j,:) - 4.0d0*var(i-1,j,:) + var(i-2,j,:))
      dvardeta = 0.5d0*(var(i,j+1,:) - var(i,j-1,:))
      gradVar(i,j,1:4) = 1.0d0/grid%elem(i,j)%detJ*(dvardxi*grid%elem(i,j)%dydeta-dvardeta*grid%elem(i,j)%dydxi)
      gradVar(i,j,5:8) = 1.0d0/grid%elem(i,j)%detJ*(dvardeta*grid%elem(i,j)%dxdxi-dvardxi*grid%elem(i,j)%dxdeta)

    end do

  end subroutine compute_gradient

  !-----------------------------------------------------------------------------
  ! Subroutine for computing the gradient velocity and temperature for interior
  ! faces using the Green-Gauss approach.
  !
  ! Process is relatively straightforward.  Need to create a new control
  ! volume whose center is the midpoint at which the gradient is to be computed.
  ! After that, need to find the value of the variables at the faces.  Two of
  ! the midpoints on the faces will be at cell centers, so those values can be
  ! taken directly.  The other two will be located at nodes and will need to
  ! be computed by averaging the values of the variables from the four
  ! adjacent cells.  After that, the gradient at the face can be computed as
  ! grad(U) = 1/Omega' ( Sum_i U_i dS_i ) where dS_i is the length of the face
  ! multiplied against the face normal.
  !
  ! Important side note:  The viscosity must be found by averaging the values 
  ! at the two cell centers.
  !-----------------------------------------------------------------------------
  subroutine compute_face_gradients(grid)
    implicit none
    type(mesh), intent(in) :: grid

  end subroutine compute_face_gradients

end module grad
