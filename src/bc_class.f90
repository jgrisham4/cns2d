!==========================================================
! This module contains a class for boundary conditions
!==========================================================

module bc_class
  use mesh_class
  implicit none
  private
  public :: bc,apply

  !--------------------------------------------------------
  ! bc class definition
  !--------------------------------------------------------
  type bc
    character (len=30)     :: bcname   ! Name of boundary condition
    integer                :: bindex   ! Index for boundary
  end type bc

  contains

    !------------------------------------------------------
    ! Subroutine for applying the boundary conditions
    ! to the solution
    !------------------------------------------------------
    subroutine apply(this,grid,bidx)
      implicit none
      type(bc), intent(inout)   :: this
      type(mesh), intent(inout) :: grid
      integer                   :: i,j

      ! Checking names of boundary conditions
      if (this%bcname.eq."periodic") then

        ! Copying state
        if (this%bindex.eq.1) then

          do i=1,grid%nelemi

            ! Copying data to ghost cells along the bottom
            grid%elem(i, 0)%u = grid%elem(i,grid%elemj)%u
            grid%elem(i,-1)%u = grid%elem(i,grid%elemj-1)%u
            grid%elem(i, 0)%w = grid%elem(i,grid%elemj)%w
            grid%elem(i,-1)%w = grid%elem(i,grid%elemj-1)%w
            grid%elem(i, 0)%w = grid%elem(i,grid%elemj)%w
            grid%elem(i,-1)%w = grid%elem(i,grid%elemj-1)%w

          end do


        else if (this%bindex.eq.2) then

        else if (this%bindex.eq.3) then

        else if (this%bindex.eq.4) then

        else

        end if

      else if (this%bcname.eq."farfield") then

      else if (this%bcname.eq."slip wall") then

      end if

    end subroutine apply

end module bc_class
