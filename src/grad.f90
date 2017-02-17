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
  use utils,      only : nvec
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
  !
  ! THE BOUNDARY CONDITIONS MUST BE ENFORCED PRIOR TO CALLING THIS SUBROUTINE
  ! OTHERWISE, THE GRADIENTS ALONG THE BOUNDARIES WILL BE WRONG.
  !-----------------------------------------------------------------------------
  subroutine compute_face_gradients(grid)
    implicit none
    type(mesh), intent(in) :: grid
    double precision :: dS(4),n(2,4),x(4),y(4),xm(4),ym(4)
    double precision :: f(3,4),wtmp(4),T

    ! Blazek says that the values in the corner ghost cells should be
    ! set using averages of the corresponding cells above or below and
    ! to the right or left.  This should work for all boundary conditions
    ! except wall or symmetry, in which case, the boundary should be
    ! effectively extended into the ghost cells.  As such, I believe
    ! this should be done in the bcs.f90 file.

    ! I'm going to compute all the geometric quantities from scratch
    ! to make sure I'm not making an indexing error with vertical
    ! and horizontal faces.

    ! Computing the gradients along vertical faces
    do j=1,this%grid%nelemj
      do i=1,this%grid%imax

        ! Computing node coordinates of the new element
        x(1) = 0.5d0*(this%grid%x(i-1,j)+this%grid%x(i,j))
        y(1) = 0.5d0*(this%grid%y(i-1,j)+this%grid%y(i,j))
        x(2) = 0.5d0*(this%grid%x(i,j)+this%grid%x(i+1,j))
        y(2) = 0.5d0*(this%grid%y(i,j)+this%grid%y(i+1,j))
        x(3) = 0.5d0*(this%grid%x(i,j+1)+this%grid%x(i+1,j+1))
        y(3) = 0.5d0*(this%grid%y(i,j+1)+this%grid%y(i+1,j+1))
        x(4) = 0.5d0*(this%grid%x(i,j+1)+this%grid%x(i-1,j+1))
        y(4) = 0.5d0*(this%grid%y(i,j+1)+this%grid%y(i-1,j+1))

        ! Computing the midpoints of the edges
        !xm(1) =
        ! NOT SURE THIS IS NECESSARY

        ! Computing the lengths of the edges
        dS(1) = sqrt((x(2)-x(1))**2+(y(2)-y(1))**2)
        dS(2) = sqrt((x(3)-x(2))**2+(y(3)-y(2))**2)
        dS(3) = sqrt((x(4)-x(3))**2+(y(4)-y(3))**2)
        dS(4) = sqrt((x(1)-x(4))**2+(y(1)-y(4))**2)

        ! Finding face normals
        n(:,1) = nvec(y(2)-y(1),x(1)-x(2))
        n(:,2) = nvec(y(3)-y(2),x(2)-x(3))
        n(:,3) = nvec(y(4)-y(3),x(3)-x(4))
        n(:,4) = nvec(y(1)-y(4),x(4)-x(1))

        ! Finding values at the midpoints of faces via averaging
        ! Only need velocities and temperature
        wtmp = u_to_w(this%grid%elem(i-1,j)%u,this%g)
        f(1,1) = wtmp(2)
        f(2,1) = wtmp(3)
        f(3,1) = wtmp(4)/(wtmp(1)*this%R)


        ! Computing the face-centered gradient

      end do
    end do

    ! Computing gradients along horizontal faces
    do j=1,this%grid%imax
      do i=1,this%grid%nelemi

        ! Computing the lengths of the edges

        ! Finding face normals

        ! Computing the face-centered gradient

      end do
    end do

    ! Computing gradients for boundary faces
    ! Requires construction of midpoints of edges


  end subroutine compute_face_gradients

end module grad
