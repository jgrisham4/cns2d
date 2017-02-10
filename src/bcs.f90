!===============================================================================
! This module contains a subroutine for applying boundary conditions to the
! boundary fluxes.
!
! Boundary conditions:
! 1000 - Farfield
! 1001 - Extrapolate
! 1002 - Weakly enforced slip wall
! 1003 - Strongly enforced slip wall
! 1004 - No-slip wall (only implemented for boundary 1)
! 2000 - Method of manufactured solution
!
! NOTES:
! - The no-slip wall BC is only implemented for the bottom boundary.  Any of
!   the boundary that is to the left of x=0 is set as an inviscid wall.  Any
!   of the boundary where x>=0 is set to viscous wall.
! - States in the first layer of ghost cells are set in a consistent manner
!   in the below.  This is necessary for Green-Gauss gradient computation.
!
! Author: James Grisham
! Date: 01/13/2017
! Modified: 02/09/2017
!===============================================================================

module bcs
  use solver_class, only : solver
  use mms,          only : s_continuity,s_xmom,s_ymom,s_energy,rho_e,u_e,v_e,et_e,dudx_e,dudy_e,dvdx_e,dvdy_e,dTdx_e,dTdy_e
  use flux,         only : flux_adv,flux_visc_state
  use riemann,      only : roe
  use utils,        only : u_to_w,nvec
  implicit none
  public :: apply_bcs

  contains

    !---------------------------------------------------------------------------
    ! Subroutine for applying boundary conditions
    !---------------------------------------------------------------------------
    subroutine apply_bcs(this)
      implicit none
      type(solver), intent(inout)    :: this
      double precision, dimension(4) :: wextrap,uextrap,utmp,wtmp
      double precision, dimension(4) :: duds,w1,w2,w3,w4
      double precision               :: n(4,2),length(4),area,um(4,2),dx,dy,nodes(4,2)
      double precision               :: pw,ds,rL(2),rR(2),xtmp,ytmp,gradu(2),gradv(2)
      double precision               :: u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy
      double precision               :: x1,x2,x3,x4
      double precision               :: y1,y2,y3,y4
      double precision               :: fctmp(4),fvtmp(4)
      integer                        :: i,j,k,l

      !===============================================
      ! Bottom boundary (inward pointing normal)
      !===============================================
      j = 1
      if (this%bcids(1).eq.1000) then

        ! Farfield
        do i=1,this%grid%nelemi

          ! Computing the flux
          this%grid%edges_h(i,j)%flux = -flux_adv(this%winfty,this%grid%elem(i,j)%n(:,1),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i,j-1)%u = w_to_u(this%winfty,this%g)

        end do

      else if (this%bcids(1).eq.1001) then

        ! Extrapolate
        do i=1,this%grid%nelemi

          ! Extrapolating the state to compute the flux
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_h(i,j)%flux = -flux_adv(wextrap,this%grid%elem(i,j)%n(:,1),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i,j-1)%u = uextrap

        end do

      else if (this%bcids(1).eq.1002) then

        ! Weak enforcement of slip wall
        ! This is accomplished by setting the fluxes using
        ! extrapolated pressure at the wall.
        do i=1,this%grid%nelemi

          ! Computing the flux
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i,j+1)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_h(i,j)%flux(1) = 0.0d0
          this%grid%edges_h(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,1)
          this%grid%edges_h(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,1)
          this%grid%edges_h(i,j)%flux(4) = 0.0d0

          ! Setting state in the ghost cell
          this%grid%elem(i,j-1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j-1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j-1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j-1)%u(4) =  this%grid%elem(i,j)%u(4)

        end do

      else if (this%bcids(1).eq.1003) then

        ! Strong enforcement of slip wall
        ! This is accomplished by setting the state at the ghost cells
        ! Reconstructing the interface states and solving the Riemann
        ! problem at the interface
        do i=1,this%grid%nelemi

          ! Setting state in the ghost cells
          this%grid%elem(i,j-1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j-1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j-1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j-1)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i,j-2)%u(1) =  this%grid%elem(i,j+1)%u(1)
          this%grid%elem(i,j-2)%u(2) = -this%grid%elem(i,j+1)%u(2)
          this%grid%elem(i,j-2)%u(3) = -this%grid%elem(i,j+1)%u(3)
          this%grid%elem(i,j-2)%u(4) =  this%grid%elem(i,j+1)%u(4)

          ! Finding the gradient using first-order accurate differences
          ds = sqrt((this%grid%elem(i,j-1)%xc - this%grid%elem(i,j-2)%xc)**2 + &
                    (this%grid%elem(i,j-1)%yc - this%grid%elem(i,j-2)%yc)**2)
          duds = (this%grid%elem(i,j-1)%u - this%grid%elem(i,j-2)%u)/ds

          ! Computing the left state of the interface
          rL(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j-1)%xc
          rL(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j-1)%yc
          this%grid%edges_h(i,j)%uL = this%grid%elem(i,j-1)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

          ! Computing the right state of the interface
          rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
          rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_h(i,j)%uR = this%grid%elem(i,j)%u + &
                                      this%grid%elem(i,j)%dudx*rR(1) + &
                                      this%grid%elem(i,j)%dudy*rR(2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL, &
            this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3),this%winfty)

        end do

      else if (this%bcids(1).eq.1004) then

        ! Enforcing no-slip wall boundary condition
        do i=1,this%grid%nelemi

          ! Seting state in ghost cells
          this%grid%elem(i,j-1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j-1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j-1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j-1)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i,j-2)%u(1) =  this%grid%elem(i,j+1)%u(1)
          this%grid%elem(i,j-2)%u(2) = -this%grid%elem(i,j+1)%u(2)
          this%grid%elem(i,j-2)%u(3) = -this%grid%elem(i,j+1)%u(3)
          this%grid%elem(i,j-2)%u(4) =  this%grid%elem(i,j+1)%u(4)

          ! Switching between inviscid and viscous wall for x < 0 and x >= 0, respectively
          if (this%grid%elem(i,j)%xc.lt.0.0d0) then

            ! Finding the gradient using first-order accurate differences
            ds = sqrt((this%grid%elem(i,j-1)%xc - this%grid%elem(i,j-2)%xc)**2 + &
                      (this%grid%elem(i,j-1)%yc - this%grid%elem(i,j-2)%yc)**2)
            duds = (this%grid%elem(i,j-1)%u - this%grid%elem(i,j-2)%u)/ds

            ! Computing the left state of the interface
            rL(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j-1)%xc
            rL(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j-1)%yc
            this%grid%edges_h(i,j)%uL = this%grid%elem(i,j-1)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

            ! Computing the right state of the interface
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc
            this%grid%edges_h(i,j)%uR = this%grid%elem(i,j)%u + &
                                        this%grid%elem(i,j)%dudx*rR(1) + &
                                        this%grid%elem(i,j)%dudy*rR(2)

            ! Solving the Riemann problem at the interface
            this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL, &
              this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3),this%winfty)

          else

            ! Finding the gradient using first-order accurate differences
            ds = sqrt((this%grid%elem(i,j-1)%xc - this%grid%elem(i,j-2)%xc)**2 + &
                      (this%grid%elem(i,j-1)%yc - this%grid%elem(i,j-2)%yc)**2)
            duds = (this%grid%elem(i,j-1)%u - this%grid%elem(i,j-2)%u)/ds

            ! Computing the left state of the interface
            rL(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j-1)%xc
            rL(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j-1)%yc
            this%grid%edges_h(i,j)%uL = this%grid%elem(i,j-1)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

            ! Computing the right state of the interface
            rR(1) = this%grid%edges_h(i,j)%xm - this%grid%elem(i,j)%xc
            rR(2) = this%grid%edges_h(i,j)%ym - this%grid%elem(i,j)%yc
            this%grid%edges_h(i,j)%uR = this%grid%elem(i,j)%u + &
                                        this%grid%elem(i,j)%dudx*rR(1) + &
                                        this%grid%elem(i,j)%dudy*rR(2)

            ! Solving the Riemann problem at the interface
            this%grid%edges_h(i,j)%flux = roe(this%grid%edges_h(i,j)%uL, &
            this%grid%edges_h(i,j)%uR,this%grid%elem(i,j-1)%n(:,3),this%winfty)

            ! Computing the gradients in the ghost cell using Green's theorem
            if (i.ne.1) then

              ! Trailing edge (assuming extrapolate conditions are used)
              if (i.eq.this%grid%nelemi) then
                this%grid%elem(i+1,j)%u   = this%grid%elem(i,j)%u
                this%grid%elem(i+1,j-1)%u = this%grid%elem(i,j-1)%u
                this%grid%elem(i+1,j-2)%u = this%grid%elem(i,j-2)%u
              end if

              ! Figuring out the nodes of the auxiliary volume
              x1 = 0.5d0*(this%grid%x(i,j)   + this%grid%x(i,j-1))
              y1 = 0.5d0*(this%grid%y(i,j)   + this%grid%y(i,j-1))
              x2 = 0.5d0*(this%grid%x(i+1,j) + this%grid%x(i+1,j-1))
              y2 = 0.5d0*(this%grid%y(i+1,j) + this%grid%y(i+1,j-1))
              x3 = 0.5d0*(this%grid%x(i+1,j) + this%grid%x(i+1,j+1))
              y3 = 0.5d0*(this%grid%y(i+1,j) + this%grid%y(i+1,j+1))
              x4 = 0.5d0*(this%grid%x(i,j)   + this%grid%x(i,j+1))
              y4 = 0.5d0*(this%grid%y(i,j)   + this%grid%y(i,j+1))

              ! Computing the face normals of the auxiliary control volume
              n(1,:) = nvec(y2-y1,x1-x2)
              n(2,:) = nvec(y3-y2,x2-x3)
              n(3,:) = nvec(y4-y3,x3-x4)
              n(4,:) = nvec(y1-y4,x4-x1)

              ! Computing the length of each face
              length(1) = sqrt((x2-x1)**2+(y2-y1)**2)
              length(2) = sqrt((x3-x2)**2+(y3-y2)**2)
              length(3) = sqrt((x4-x3)**2+(y4-y3)**2)
              length(4) = sqrt((x1-x4)**2+(y1-y4)**2)

              ! Computing the area of the auxiliary control volume
              area = 0.5d0*((x1-x3)*(y2-y4)+(x4-x2)*(y1-y3))

              ! Velocities at the bottom midpoint (face 1)
              um(1,1)  = this%grid%elem(i,j-1)%u(2)/this%grid%elem(i,j-1)%u(1)
              um(1,2)  = this%grid%elem(i,j-1)%u(3)/this%grid%elem(i,j-1)%u(1)

              ! Velocities at the right midpoint (face 2)
              w1 = u_to_w(this%grid%elem(i,j-1)%u,this%g)
              w2 = u_to_w(this%grid%elem(i+1,j-1)%u,this%g)
              w3 = u_to_w(this%grid%elem(i+1,j)%u,this%g)
              w4 = u_to_w(this%grid%elem(i,j)%u,this%g)
              um(2,1) = 0.25d0*(w1(2) + w2(2) + w3(2) + w4(2))
              um(2,2) = 0.25d0*(w1(3) + w2(3) + w3(3) + w4(3))

              ! Velocities at the top midpoint (face 3)
              um(3,1) = this%grid%elem(i,j)%u(2)/this%grid%elem(i,j)%u(1)
              um(3,2) = this%grid%elem(i,j)%u(3)/this%grid%elem(i,j)%u(1)

              ! Velocities at the left midpoint (face 4)
              w2 = w1
              w3 = w4
              w1 = u_to_w(this%grid%elem(i-1,j-1)%u,this%g)
              w4 = u_to_w(this%grid%elem(i-1,j)%u,this%g)
              um(4,1) = 0.25d0*(w1(2) + w2(2) + w3(2) + w4(2))
              um(4,2) = 0.25d0*(w1(3) + w2(3) + w3(3) + w4(3))

              ! Computing the gradient of velocity
              gradu(:) = 0.0d0
              gradv(:) = 0.0d0
              do k=1,4
                gradu = gradu + 1.0d0/area*(um(k,1)*n(k,:)*length(k))
                gradv = gradv + 1.0d0/area*(um(k,2)*n(k,:)*length(k))
              end do

              ! Finding average T at the face
              T = w3(4)/(w3(1)*this%R)

              ! Now have the velocity gradient at the boundary
              ! Can now compute the viscous fluxes
              this%grid%edges_h(i,j)%flux = this%grid%edges_h(i,j)%flux - &
                flux_visc_state(0.0d0,0.0d0,T,gradu(1),gradu(2),gradv(1),gradv(2),0.0d0,0.0d0,-this%grid%elem(i,j)%n(:,1),this%g,this%R)

            else

              ! Near the leading edge
              print *, "Warning: No-slip BC is not enforced here..."

            end if

          end if
        end do

      else if (this%bcids(1).eq.2000) then

        ! Method of manufactured solution BC
        j = 1
        do i=1,this%grid%nelemi

          ! Computing the exact solution at the midpoint of face
          xtmp = this%grid%edges_h(i,j)%xm
          ytmp = this%grid%edges_h(i,j)%ym
          u    = u_e(xtmp,ytmp)
          v    = v_e(xtmp,ytmp)
          T    = (this%g-1.0d0)/this%R*(et_e(xtmp,ytmp)-0.5d0*(u**2+v**2))
          dudx = dudx_e(xtmp,ytmp)
          dudy = dudy_e(xtmp,ytmp)
          dvdx = dvdx_e(xtmp,ytmp)
          dvdy = dvdy_e(xtmp,ytmp)
          dTdx = dTdx_e(xtmp,ytmp)
          dTdy = dTdy_e(xtmp,ytmp)

          ! Setting the conserved variables at the face
          utmp(1) = rho_e(xtmp,ytmp)
          utmp(2) = rho_e(xtmp,ytmp)*u
          utmp(3) = rho_e(xtmp,ytmp)*v
          utmp(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

          ! Computing the primitive variables at the face and computing the flux
          wtmp = u_to_w(utmp,this%g)
          this%grid%edges_h(i,j)%flux = flux_adv(wtmp,-this%grid%elem(i,j)%n(:,1),this%g) - &
            flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,-this%grid%elem(i,j)%n(:,1),this%g,this%R)

          ! Setting state in the ghost cells
          xtmp = this%grid%elem(i,j-1)%xc
          ytmp = this%grid%elem(i,j-1)%yc
          this%grid%elem(i,j-1)%u(1) = rho_e(xtmp,ytmp)
          this%grid%elem(i,j-1)%u(2) = rho_e(xtmp,ytmp)*u_e(xtmp,ytmp)
          this%grid%elem(i,j-1)%u(3) = rho_e(xtmp,ytmp)*v_e(xtmp,ytmp)
          this%grid%elem(i,j-1)%u(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

        end do

      else

        print *, "Boundary condition ", this%bcids(1), " for bottom not recognized."
        stop

      end if

      !===============================================
      ! Right boundary
      !===============================================
      i = this%grid%nelemi
      if (this%bcids(2).eq.1000) then

        ! Farfield
        do j=1,this%grid%nelemj

          ! Computing the flux
          this%grid%edges_v(i+1,j)%flux = flux_adv(this%winfty,this%grid%elem(i,j)%n(:,2),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i+1,j)%u = w_to_u(this%grid%elem(i,j)%u,this%g)

        end do

      else if (this%bcids(2).eq.1001) then

        ! Extrapolate
        do j=1,this%grid%nelemj

          ! Computing the flux
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_v(i+1,j)%flux = flux_adv(wextrap,this%grid%elem(i,j)%n(:,2),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i+1,j)%u = uextrap

        end do

      else if (this%bcids(2).eq.1002) then

        ! Slip wall
        do j=1,this%grid%nelemj

          ! Computing the flux
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i-1,j)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_v(i+1,j)%flux(1) = 0.0d0
          this%grid%edges_v(i+1,j)%flux(2) = pw*this%grid%elem(i,j)%n(1,2)
          this%grid%edges_v(i+1,j)%flux(3) = pw*this%grid%elem(i,j)%n(2,2)
          this%grid%edges_v(i+1,j)%flux(4) = 0.0d0

          ! Setting the state in the ghost cells
          this%grid%elem(i+1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i+1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i+1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i+1,j)%u(4) =  this%grid%elem(i,j)%u(4)

        end do

      else if (this%bcids(2).eq.1003) then

        ! Slip wall strongly enforced
        do j=1,this%grid%nelemj

          ! Setting state in ghost cells
          this%grid%elem(i+1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i+1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i+1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i+1,j)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i+2,j)%u(1) =  this%grid%elem(i-1,j)%u(1)
          this%grid%elem(i+2,j)%u(2) = -this%grid%elem(i-1,j)%u(2)
          this%grid%elem(i+2,j)%u(3) = -this%grid%elem(i-1,j)%u(3)
          this%grid%elem(i+2,j)%u(4) =  this%grid%elem(i-1,j)%u(4)

          ! Computing the gradient in the ghost cells using
          ! first-order accurate differences
          ds = sqrt((this%grid%edges_v(i+1,j)%xm - this%grid%elem(i+1,j)%xc)**2 + &
                    (this%grid%edges_v(i+1,j)%ym - this%grid%elem(i+1,j)%yc)**2)
          duds = (this%grid%elem(i+1,j)%u - this%grid%elem(i+2,j)%u)/ds

          ! Reconstructing the state on the interior
          rL(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i,j)%xc
          rL(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_v(i+1,j)%uL = this%grid%elem(i,j)%u + this%grid%elem(i,j)%dudx*rL(1) + &
            this%grid%elem(i,j)%dudy*rL(2)

          ! Reconstructing the state on the exterior
          rR(1) = this%grid%edges_v(i+1,j)%xm - this%grid%elem(i+1,j)%xc
          rR(2) = this%grid%edges_v(i+1,j)%ym - this%grid%elem(i+1,j)%yc
          this%grid%edges_v(i+1,j)%uR = this%grid%elem(i+1,j)%u + duds*sqrt(rR(1)**2+rR(2)**2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_v(i+1,j)%flux = roe(this%grid%edges_v(i+1,j)%uL, &
            this%grid%edges_v(i+1,j)%uR,this%grid%elem(i,j)%n(:,2),this%winfty)

        end do

      else if (this%bcids(2).eq.2000) then

        ! Method of manufactured solution BC
        do j=1,this%grid%nelemj

          ! Computing exact solution at midpoint of face
          xtmp = this%grid%edges_v(i+1,j)%xm
          ytmp = this%grid%edges_v(i+1,j)%ym
          u    = u_e(xtmp,ytmp)
          v    = v_e(xtmp,ytmp)
          T    = (this%g-1.0d0)/this%R*(et_e(xtmp,ytmp)-0.5d0*(u**2+v**2))
          dudx = dudx_e(xtmp,ytmp)
          dudy = dudy_e(xtmp,ytmp)
          dvdx = dvdx_e(xtmp,ytmp)
          dvdy = dvdy_e(xtmp,ytmp)
          dTdx = dTdx_e(xtmp,ytmp)
          dTdy = dTdy_e(xtmp,ytmp)

          ! Setting conserved variables at the face
          utmp(1) = rho_e(xtmp,ytmp)
          utmp(2) = rho_e(xtmp,ytmp)*u
          utmp(3) = rho_e(xtmp,ytmp)*v
          utmp(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

          ! Computing primitive variables at the face and computing the flux
          wtmp = u_to_w(utmp,this%g)
          this%grid%edges_v(i+1,j)%flux = flux_adv(wtmp,this%grid%elem(i,j)%n(:,2),this%g) - &
            flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,this%grid%elem(i,j)%n(:,2),this%g,this%R)

          ! Setting state in the ghost cell
          xtmp = this%grid%elem(i+1,j)%xc
          ytmp = this%grid%elem(i+1,j)%yc
          this%grid%elem(i+1,j)%u(1) = rho_e(xtmp,ytmp)
          this%grid%elem(i+1,j)%u(2) = rho_e(xtmp,ytmp)*u_e(xtmp,ytmp)
          this%grid%elem(i+1,j)%u(3) = rho_e(xtmp,ytmp)*v_e(xtmp,ytmp)
          this%grid%elem(i+1,j)%u(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

        end do

      else

        print *, "Boundary condition ", this%bcids(2), " for right not recognized."
        stop

      end if

      !===============================================
      ! Top boundary
      !===============================================
      j = this%grid%nelemj
      if (this%bcids(3).eq.1000) then

        ! Farfield
        do i=1,this%grid%nelemi

          ! Computing flux
          this%grid%edges_h(i,j+1)%flux = flux_adv(this%winfty,this%grid%elem(i,j)%n(:,3),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i,j+1)%u = w_to_u(this%winfty,this%g)

        end do

      else if (this%bcids(3).eq.1001) then

        ! Extrapolate
        do i=1,this%grid%nelemi

          ! Computing the flux
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_h(i,j+1)%flux = flux_adv(wextrap,this%grid%elem(i,j)%n(:,3),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i,j+1)%u = uextrap

        end do

      else if (this%bcids(3).eq.1002) then

        ! Weak enforcement of slip wall
        do i=1,this%grid%nelemi

          ! Computing the flux
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i,j-1)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_h(i,j+1)%flux(1) = 0.0d0
          this%grid%edges_h(i,j+1)%flux(2) = pw*this%grid%elem(i,j)%n(1,3)
          this%grid%edges_h(i,j+1)%flux(3) = pw*this%grid%elem(i,j)%n(2,3)
          this%grid%edges_h(i,j+1)%flux(4) = 0.0d0

          ! Setting the state in the ghost cell
          this%grid%elem(i,j+1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j+1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j+1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j+1)%u(4) =  this%grid%elem(i,j)%u(4)

        end do

      else if (this%bcids(3).eq.1003) then

        ! Strong enforcement of slip wall
        ! This is accomplished by setting the state at the ghost cells
        ! Reconstructing the interface states and solving the Riemann
        ! problem at the interface
        do i=1,this%grid%nelemi

          ! Setting state in the ghost cells
          this%grid%elem(i,j+1)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i,j+1)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i,j+1)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i,j+1)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i,j+2)%u(1) =  this%grid%elem(i,j-1)%u(1)
          this%grid%elem(i,j+2)%u(2) = -this%grid%elem(i,j-1)%u(2)
          this%grid%elem(i,j+2)%u(3) = -this%grid%elem(i,j-1)%u(3)
          this%grid%elem(i,j+2)%u(4) =  this%grid%elem(i,j-1)%u(4)

          ! Finding the gradient using first-order accurate differences
          ds = sqrt((this%grid%elem(i,j+1)%xc - this%grid%elem(i,j+2)%xc)**2 + &
                    (this%grid%elem(i,j+1)%yc - this%grid%elem(i,j+2)%yc)**2)
          duds = (this%grid%elem(i,j+1)%u - this%grid%elem(i,j+2)%u)/ds

          ! o------------o
          ! |            | j+2 -- ghost
          ! |            |
          ! o------------o
          ! |            | j+1 -- ghost
          ! | R          |
          ! o------------o
          ! | L          | j
          ! |            |
          ! o------------o
          ! Computing the left state of the interface
          rL(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j)%xc
          rL(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_h(i,j+1)%uL = this%grid%elem(i,j)%u + &
                                        this%grid%elem(i,j)%dudx*rL(1) + &
                                        this%grid%elem(i,j)%dudy*rL(2)

          ! Computing the right state of the interface
          rR(1) = this%grid%edges_h(i,j+1)%xm - this%grid%elem(i,j+1)%xc
          rR(2) = this%grid%edges_h(i,j+1)%ym - this%grid%elem(i,j+1)%yc
          this%grid%edges_h(i,j+1)%uR = this%grid%elem(i,j+1)%u + duds*sqrt(rR(1)**2 + rR(2)**2)

          ! Solving the Riemann problem at the interface
          this%grid%edges_h(i,j+1)%flux = roe(this%grid%edges_h(i,j+1)%uL, &
            this%grid%edges_h(i,j+1)%uR,this%grid%elem(i,j)%n(:,3),this%winfty)

        end do

      else if (this%bcids(3).eq.2000) then

        ! Method of manufactured solution BC
        do i=1,this%grid%nelemi

          ! Computing the exact solution at the midpoint of the face
          xtmp = this%grid%edges_h(i,j+1)%xm
          ytmp = this%grid%edges_h(i,j+1)%ym
          u    = u_e(xtmp,ytmp)
          v    = v_e(xtmp,ytmp)
          T    = (this%g-1.0d0)/this%R*(et_e(xtmp,ytmp)-0.5d0*(u**2+v**2))
          dudx = dudx_e(xtmp,ytmp)
          dudy = dudy_e(xtmp,ytmp)
          dvdx = dvdx_e(xtmp,ytmp)
          dvdy = dvdy_e(xtmp,ytmp)
          dTdx = dTdx_e(xtmp,ytmp)
          dTdy = dTdy_e(xtmp,ytmp)

          ! Setting the conserved variables at the face
          utmp(1) = rho_e(xtmp,ytmp)
          utmp(2) = rho_e(xtmp,ytmp)*u
          utmp(3) = rho_e(xtmp,ytmp)*v
          utmp(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

          ! Computing primitive variables at the face and computing the flux
          wtmp = u_to_w(utmp,this%g)
          this%grid%edges_h(i,j+1)%flux = flux_adv(wtmp,this%grid%elem(i,j)%n(:,3),this%g) - &
            flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,this%grid%elem(i,j)%n(:,3),this%g,this%R)

          ! Setting the state in the ghost cell
          xtmp = this%grid%elem(i,j+1)%xc
          ytmp = this%grid%elem(i,j+1)%yc
          this%grid%elem(i,j+1)%u(1) = rho_e(xtmp,ytmp)
          this%grid%elem(i,j+1)%u(2) = rho_e(xtmp,ytmp)*u_e(xtmp,ytmp)
          this%grid%elem(i,j+1)%u(3) = rho_e(xtmp,ytmp)*v_e(xtmp,ytmp)
          this%grid%elem(i,j+1)%u(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

        end do

      else

        print *, "Boundary condition ", this%bcids(3), " for top not recognized."
        stop

      end if

      !===============================================
      ! Left boundary (inward pointing normal)
      !===============================================
      i = 1
      if (this%bcids(4).eq.1000) then

        ! Farfield
        do j=1,this%grid%nelemj

          ! Computing the flux
          this%grid%edges_v(i,j)%flux = -flux_adv(this%winfty,this%grid%elem(i,j)%n(:,4),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i-1,j)%u = w_to_u(this%winfty,this%g)

        end do

      else if (this%bcids(4).eq.1001) then

        ! Extrapolate
        do j=1,this%grid%nelemj

          ! Computing the flux
          uextrap = this%grid%elem(i,j)%u
          wextrap = u_to_w(uextrap,this%g)
          this%grid%edges_v(i,j)%flux = -flux_adv(wextrap,this%grid%elem(i,j)%n(:,4),this%g)

          ! Setting the state in the ghost cell
          this%grid%elem(i-1,j)%u = uextrap

        end do

      else if (this%bcids(4).eq.1002) then

        ! Weakly enforced slip wall
        do j=1,this%grid%nelemj

          ! Computing the flux
          wtmp = u_to_w(this%grid%elem(i,j)%u,this%g)
          pw = 3.0d0/2.0d0*wtmp(4)
          wtmp = u_to_w(this%grid%elem(i+1,j)%u,this%g)
          pw = pw - 0.5d0*wtmp(4)
          this%grid%edges_v(i,j)%flux(1) = 0.0d0
          this%grid%edges_v(i,j)%flux(2) = -pw*this%grid%elem(i,j)%n(1,4)
          this%grid%edges_v(i,j)%flux(3) = -pw*this%grid%elem(i,j)%n(2,4)
          this%grid%edges_v(i,j)%flux(4) = 0.0d0

          ! Setting the state in the ghost cell
          this%grid%elem(i-1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i-1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i-1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i-1,j)%u(4) =  this%grid%elem(i,j)%u(4)

        end do

      else if (this%bcids(4).eq.1003) then

        ! Strongly enforced slip wall
        do j=1,this%grid%nelemj

          ! Assigning states in the ghost cells
          this%grid%elem(i-1,j)%u(1) =  this%grid%elem(i,j)%u(1)
          this%grid%elem(i-1,j)%u(2) = -this%grid%elem(i,j)%u(2)
          this%grid%elem(i-1,j)%u(3) = -this%grid%elem(i,j)%u(3)
          this%grid%elem(i-1,j)%u(4) =  this%grid%elem(i,j)%u(4)
          this%grid%elem(i-2,j)%u(1) =  this%grid%elem(i+1,j)%u(1)
          this%grid%elem(i-2,j)%u(2) = -this%grid%elem(i+1,j)%u(2)
          this%grid%elem(i-2,j)%u(3) = -this%grid%elem(i+1,j)%u(3)
          this%grid%elem(i-2,j)%u(4) =  this%grid%elem(i+1,j)%u(4)

          ! Computing gradient in the ghost cell
          ds = sqrt((this%grid%elem(i-1,j)%xc - this%grid%elem(i-2,j)%xc)**2 + &
                    (this%grid%elem(i-1,j)%yc - this%grid%elem(i-2,j)%yc)**2)
          duds = (this%grid%elem(i-1,j)%u - this%grid%elem(i-2,j)%u)/ds

          ! Reconstructing the exterior state
          rL(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i-1,j)%xc
          rL(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i-1,j)%yc
          this%grid%edges_v(i,j)%uL = this%grid%elem(i-1,j)%u + duds*sqrt(rL(1)**2 + rL(2)**2)

          ! Reconstructing the interior state
          rR(1) = this%grid%edges_v(i,j)%xm - this%grid%elem(i,j)%xc
          rR(2) = this%grid%edges_v(i,j)%ym - this%grid%elem(i,j)%yc
          this%grid%edges_v(i,j)%uR = this%grid%elem(i,j)%u + &
                                      this%grid%elem(i,j)%dudx*rR(1) + &
                                      this%grid%elem(i,j)%dudy*rR(2)

          ! Solving the Riemann problem
          this%grid%edges_v(i,j)%flux = roe(this%grid%edges_v(i,j)%uL, &
            this%grid%edges_v(i,j)%uR,this%grid%elem(i-1,j)%n(:,2),this%winfty)

        end do

      else if (this%bcids(4).eq.2000) then

        ! Method of manufactured solution BC
        do j=1,this%grid%nelemj

          ! Computing the exact solution at the midpoint of face
          xtmp = this%grid%edges_v(i,j)%xm
          ytmp = this%grid%edges_v(i,j)%ym
          u    = u_e(xtmp,ytmp)
          v    = v_e(xtmp,ytmp)
          T    = (this%g-1.0d0)/this%R*(et_e(xtmp,ytmp)-0.5d0*(u**2+v**2))
          dudx = dudx_e(xtmp,ytmp)
          dudy = dudy_e(xtmp,ytmp)
          dvdx = dvdx_e(xtmp,ytmp)
          dvdy = dvdy_e(xtmp,ytmp)
          dTdx = dTdx_e(xtmp,ytmp)
          dTdy = dTdy_e(xtmp,ytmp)

          ! Setting the conserved variables at the face
          utmp(1) = rho_e(xtmp,ytmp)
          utmp(2) = rho_e(xtmp,ytmp)*u
          utmp(3) = rho_e(xtmp,ytmp)*v
          utmp(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

          ! Computing the primitive variables at the face and computing the flux
          wtmp = u_to_w(utmp,this%g)
          this%grid%edges_v(i,j)%flux = flux_adv(wtmp,-this%grid%elem(i,j)%n(:,4),this%g) - &
            flux_visc_state(u,v,T,dudx,dudy,dvdx,dvdy,dTdx,dTdy,-this%grid%elem(i,j)%n(:,4),this%g,this%R)

          ! Setting the state in the ghost cell
          xtmp = this%grid%elem(i-1,j)%xc
          ytmp = this%grid%elem(i-1,j)%yc
          this%grid%elem(i-1,j)%u(1) = rho_e(xtmp,ytmp)
          this%grid%elem(i-1,j)%u(2) = rho_e(xtmp,ytmp)*u_e(xtmp,ytmp)
          this%grid%elem(i-1,j)%u(3) = rho_e(xtmp,ytmp)*v_e(xtmp,ytmp)
          this%grid%elem(i-1,j)%u(4) = rho_e(xtmp,ytmp)*et_e(xtmp,ytmp)

        end do

      else

        print *, "Boundary condition ", this%bcids(4), " for left not recognized."
        stop

      end if

      ! Filling in values in the ghost cells
      ! Bottom left
      i = 0
      j = 0
      if (this%bcids(1)=="1004") then
        this%grid%elem(i,j)%u(1) =  this%grid%elem(i,j+1)%u(1)
        this%grid%elem(i,j)%u(2) = -this%grid%elem(i,j+1)%u(2)
        this%grid%elem(i,j)%u(3) = -this%grid%elem(i,j+1)%u(3)
        this%grid%elem(i,j)%u(4) =  this%grid%elem(i,j+1)%u(4)
      else
        this%grid%elem(i,j)%u = 0.5d0*(this%grid%elem(i+1,j)%u + &
          this%grid%elem(i,j+1)%u)
      endif

      ! Bottom right
      i = this%grid%imax+1
      j = 0
      if (this%bcids(1)=="1004") then
        this%grid%elem(i,j)%u(1) =  this%grid%elem(i,j+1)%u(1)
        this%grid%elem(i,j)%u(2) = -this%grid%elem(i,j+1)%u(2)
        this%grid%elem(i,j)%u(3) = -this%grid%elem(i,j+1)%u(3)
        this%grid%elem(i,j)%u(4) =  this%grid%elem(i,j+1)%u(4)
      else
        this%grid%elem(i,j)%u = 0.5d0*(this%grid%elem(i-1,j)%u + &
          this%grid%elem(i,j+1)%u)
      endif

      ! Top right
      i = this%grid%imax+1
      j = this%grid%jmax+1
      this%grid%elem(i,j)%u = 0.5d0*(this%grid%elem(i-1,j)%u + &
        this%grid%elem(i,j-1)%u)

      ! Top left
      i = 0
      j = this%grid%jmax+1
      this%grid%elem(i,j)%u = 0.5d0*(this%grid%elem(i+1,j)%u + &
        this%grid%elem(i,j-1)%u)

    end subroutine apply_bcs

end module bcs
