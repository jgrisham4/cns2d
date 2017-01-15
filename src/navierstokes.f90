!===============================================================================
! This module contains the subroutines for a simple Navier-Stokes solver.
!
! Author: James Grisham
! Date: 01/13/2017
!===============================================================================

module navierstokes
  use solver_class, only : solver
  use utils,        only : u_to_w
  use flux,         only : flux_visc
  use grad,         only : compute_gradient
  use euler,        only : residual_inv,residual_inv_fo
  implicit none
  public :: residual_visc,residual_visc_fo

  contains

    !---------------------------------------------------------------------------
    ! Subroutine for computing the residual for Navier-Stokes
    ! equations
    !---------------------------------------------------------------------------
    subroutine residual_visc(this,resid)
      implicit none
      type(solver), intent(inout)     :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4),phi(4)
      double precision, dimension(4)  :: fx,fy,uextrap,wextrap,wtmp
      double precision, dimension(2)  :: rL,rR,r
      double precision, allocatable   :: u(:,:,:),gradU(:,:,:),umax(:,:,:),umin(:,:,:)
      double precision                :: src(4),xtmp,ytmp  ! only used for MMS
      integer                         :: i,j,k,l,err

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in residual_visc."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in residual_visc."
        stop
      end if

      ! Finding inviscid part of the residual
      call residual_inv(this,resid)

      ! Computing gradient of temperature and velocity
      ! dT/dx_ij = gradU(i,j,1)
      ! dT/dy_ij = gradU(i,j,5)
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%w = u_to_w(this%grid%elem(i,j)%u,this%g)
          u(i,j,1) = this%grid%elem(i,j)%w(4)/(this%R*this%grid%elem(i,j)%w(1))
          u(i,j,2) = this%grid%elem(i,j)%w(2)
          u(i,j,3) = this%grid%elem(i,j)%w(3)
          u(i,j,4) = 0.0d0
        end do
      end do
      call compute_gradient(this%grid,u,gradU)

      ! Copying gradients of temperature and velocity to element objects
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%dTdx = gradU(i,j,1)
          this%grid%elem(i,j)%dTdy = gradU(i,j,5)
          this%grid%elem(i,j)%dVdx = gradU(i,j,2:3)  ! Vector which holds {du/dx, dv/dx}
          this%grid%elem(i,j)%dVdy = gradU(i,j,6:7)  ! Vector which holds {du/dy, dv/dy}
        end do
      end do

      ! Computing viscous fluxes for vertical internal faces
      ! The flux for each edge has already been set by the residual_inv subroutine
      ! Just need to subtract the viscous flux from the inviscid flux
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi
          this%grid%edges_v(i,j)%flux = this%grid%edges_v(i,j)%flux - &
            flux_visc(this%grid%elem(i-1,j),this%grid%elem(i,j), &
                      this%grid%elem(i,j)%n(:,2),this%g,this%R)
        end do
      end do

      ! Computing viscous fluxes for the horizontal internal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j)%flux = this%grid%edges_h(i,j)%flux - &
            flux_visc(this%grid%elem(i,j-1),this%grid%elem(i,j), &
                      this%grid%elem(i,j)%n(:,3),this%g,this%R)
        end do
      end do

      ! Computing residual
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          resid(i,j,:) = -1.0d0/(this%grid%elem(i,j)%area)*                 &
            (-this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length +   &
            this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length - &
            this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length +     &
            this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length)
        end do
      end do

      ! Computing residual - METHOD OF MANUFACTURED SOLUTIONS
      !do j=1,this%grid%nelemj
      !  do i=1,this%grid%nelemi

      !    ! Computing source terms
      !    xtmp   = this%grid%elem(i,j)%xc
      !    ytmp   = this%grid%elem(i,j)%yc
      !    src(1) = s_continuity(xtmp,ytmp)
      !    src(2) = s_xmom(xtmp,ytmp)
      !    src(3) = s_ymom(xtmp,ytmp)
      !    src(4) = s_energy(xtmp,ytmp)

      !    ! Forming residual
      !    resid(i,j,:) = -1.0d0/(this%grid%elem(i,j)%area)*                 &
      !      (-this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length +   &
      !      this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length - &
      !      this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length +     &
      !      this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length - &
      !      src*this%grid%elem(i,j)%area)

      !  end do
      !end do

    end subroutine residual_visc

    !---------------------------------------------------------------------------
    ! Subroutine for computing the residual for Navier-Stokes
    ! equations using first-order spatial accuracy.
    !---------------------------------------------------------------------------
    subroutine residual_visc_fo(this,resid)
      implicit none
      type(solver), intent(inout)     :: this
      double precision, intent(inout) :: resid(:,:,:)
      double precision                :: duL(4),duR(4),phi(4)
      double precision, dimension(4)  :: fx,fy,uextrap,wextrap,wtmp
      double precision, dimension(2)  :: rL,rR,r
      double precision, allocatable   :: u(:,:,:),gradU(:,:,:),umax(:,:,:),umin(:,:,:)
      double precision                :: src(4),xtmp,ytmp  ! only used for MMS
      integer                         :: i,j,k,l,err

      ! Allocating memory for gradient of state at cell centers
      allocate(gradU(this%grid%nelemi,this%grid%nelemj,8),stat=err)
      if (err.ne.0) then
        print *, "Error: Can't allocate memory for gradU in residual_visc."
        stop
      end if
      allocate(u(this%grid%nelemi,this%grid%nelemj,4),stat=err)
      if (err.ne.0) then
        print *, "Error: can't allocate memory for u in residual_visc."
        stop
      end if

      ! Finding inviscid part of the residual
      call residual_inv_fo(this,resid)

      ! Computing gradient of temperature and velocity
      ! dT/dx_ij = gradU(i,j,1)
      ! dT/dy_ij = gradU(i,j,5)
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%w = u_to_w(this%grid%elem(i,j)%u,this%g)
          u(i,j,1) = this%grid%elem(i,j)%w(4)/(this%R*this%grid%elem(i,j)%w(1))
          u(i,j,2) = this%grid%elem(i,j)%w(2)
          u(i,j,3) = this%grid%elem(i,j)%w(3)
          u(i,j,4) = 0.0d0
        end do
      end do
      call compute_gradient(this%grid,u,gradU)

      ! Copying gradients of temperature and velocity to element objects
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%elem(i,j)%dTdx = gradU(i,j,1)
          this%grid%elem(i,j)%dTdy = gradU(i,j,5)
          this%grid%elem(i,j)%dVdx = gradU(i,j,2:3)  ! Vector which holds {du/dx, dv/dx}
          this%grid%elem(i,j)%dVdy = gradU(i,j,6:7)  ! Vector which holds {du/dy, dv/dy}
        end do
      end do

      ! Computing viscous fluxes for vertical internal faces
      ! The flux for each edge has already been set by the residual_inv subroutine
      ! Just need to subtract the viscous flux from the inviscid flux
      do j=1,this%grid%nelemj
        do i=2,this%grid%nelemi
          this%grid%edges_v(i,j)%flux = this%grid%edges_v(i,j)%flux - &
            flux_visc(this%grid%elem(i-1,j),this%grid%elem(i,j), &
                      this%grid%elem(i,j)%n(:,2),this%g,this%R)
        end do
      end do

      ! Computing viscous fluxes for the horizontal internal faces
      do j=2,this%grid%nelemj
        do i=1,this%grid%nelemi
          this%grid%edges_h(i,j)%flux = this%grid%edges_h(i,j)%flux - &
            flux_visc(this%grid%elem(i,j-1),this%grid%elem(i,j), &
                      this%grid%elem(i,j)%n(:,3),this%g,this%R)
        end do
      end do

      ! Computing residual
      do j=1,this%grid%nelemj
        do i=1,this%grid%nelemi
          resid(i,j,:) = -1.0d0/(this%grid%elem(i,j)%area)*                 &
            (-this%grid%edges_h(i,j)%flux*this%grid%edges_h(i,j)%length +   &
            this%grid%edges_h(i,j+1)%flux*this%grid%edges_h(i,j+1)%length - &
            this%grid%edges_v(i,j)%flux*this%grid%edges_v(i,j)%length +     &
            this%grid%edges_v(i+1,j)%flux*this%grid%edges_v(i+1,j)%length)
        end do
      end do

    end subroutine residual_visc_fo

end module navierstokes
