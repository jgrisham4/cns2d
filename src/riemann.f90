!*****************************************************************************
!* A collection of two-dimensional Euler numerical fluxes, Version 2 (2010),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*------------------------
!* List of Flux Functions:
!
!   Roe
!   Rotated-RHLL
!
!*------------------------
!*
!* These F90 routines were written and made available for download
!* for an educational purpose. Detailed descripstion of each numerical
!* flux can be found in the original paper, or in popular textbooks, or
!* in the (will-be-available) second volume of "I do like CFD".
!*
!* Note that all routines are not efficietly implemented for clarity; you can
!* improve the efficiecy and also you can covert it to double precision
!* version if you wish.
!*
!* NOTES: There were bugs in the Rotated-RHLL solver in version 1.
!*        It has been fixed and tested for a shock-diffraction problem.
!*
!* This file may be updated in future. (to add more flux functions.)
!*****************************************************************************

!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
module riemann
  implicit none
  private
  public :: roe, rotated_rhll
  contains
    function roe(uL, uR, nhat, winf)
      double precision :: uL(4), uR(4) ! Input: conservative variables rho*[1, u, v, E]
      double precision :: nhat(2)      ! Input: face normal vector
      double precision :: winf(4)      ! Input: Freestream primitives
      double precision :: nx, ny       ! Input: Face normal vector, [nx, ny] (Left-to-Right)
      double precision :: roe(4)       ! Output: Roe flux function (upwind)
      !Local constants
      double precision :: gam                          ! Ratio of specific heat.
      double precision :: zero, fifth, half, one, two    ! Numbers
      !Local variables
      double precision :: tx, ty                         ! Tangent vector (perpendicular to the face normal)
      double precision :: vxL, vxR, vyL, vyR             ! Velocity components.
      double precision :: rhoL, rhoR, pL, pR             ! Primitive variables.
      double precision :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
      double precision :: aL, aR, HL, HR                 ! Speeds of sound.
      double precision :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
      double precision :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
      double precision :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
      double precision :: fL(4), fR(4), diss(4)          ! Fluxes and dissipation term
      double precision, parameter :: floor_val = 0.01d0  ! Floor value used to limit nonphysical p and rho
      integer :: i, j

      !Constants.
      gam = 1.4d0
      zero = 0.0d0
      fifth = 0.2d0
      half = 0.5d0
      one = 1.0d0
      two = 2.0d0

      !Tangent vector (Do you like it? Actually, Roe flux can be implemented
      ! without any tangent vector. See "I do like CFD, VOL.1" for details.)
      nx = nhat(1)
      ny = nhat(2)
      tx = -ny
      ty = nx

      !Primitive and other variables.
      !  Left state
      rhoL = uL(1)
      vxL  = uL(2)/uL(1)
      vyL  = uL(3)/uL(1)
      vnL  = vxL*nx+vyL*ny
      vtL  = vxL*tx+vyL*ty
      if (rhoL<(floor_val*winf(1))) then
        print *, "Warning: rhoL floored in roe subroutine."
        rhoL = floor_val*winf(1)
      end if
      pL   = (gam-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )  ! Getting a negative pressure
      if (pL<(floor_val*winf(4))) then
        print *, "Warning: pL floored in roe subroutine."
        pL = floor_val*winf(4)
      end if
      aL   = sqrt(gam*pL/rhoL)
      HL   = ( uL(4) + pL ) / rhoL

      !  Right state
      rhoR = uR(1)
      vxR = uR(2)/uR(1)
      vyR = uR(3)/uR(1)
      vnR = vxR*nx+vyR*ny
      vtR = vxR*tx+vyR*ty
      if (rhoR<(floor_val*winf(1))) then
        print *, "Warning: rhoR floored in roe subroutine."
        rhoR = floor_val*winf(1)
      end if
      pR = (gam-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      if (pR<(floor_val*winf(4))) then
        print *, "Warning: pR floored in roe subroutine."
        pR = floor_val*winf(4)
      end if
      aR = sqrt(gam*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

      !First compute the Roe Averages
      RT  = sqrt(rhoR/rhoL)
      rho = RT*rhoL
      vx  = (vxL+RT*vxR)/(one+RT)
      vy  = (vyL+RT*vyR)/(one+RT)
      H   = ( HL+RT* HR)/(one+RT)
      a   = sqrt( (gam-one)*(H-half*(vx*vx+vy*vy)) )
      vn  = vx*nx+vy*ny
      vt  = vx*tx+vy*ty

      !Wave Strengths
      drho = rhoR - rhoL
      dp =   pR - pL
      dvn =  vnR - vnL
      dvt =  vtR - vtL

      dV(1) = (dp - rho*a*dvn )/(two*a*a)
      dV(2) = rho*dvt/a
      dV(3) =  drho - dp/(a*a)
      dV(4) = (dp + rho*a*dvn )/(two*a*a)

      !Wave Speed
      ws(1) = abs(vn-a)
      ws(2) = abs(vn)
      ws(3) = abs(vn)
      ws(4) = abs(vn+a)

      !Harten's Entropy Fix JCP(1983), 49, pp357-393:
      ! only for the nonlinear fields.
      dws(1) = fifth
      if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
      dws(4) = fifth
      if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

      !Right Eigenvectors
      Rv(1,1) = one
      Rv(2,1) = vx - a*nx
      Rv(3,1) = vy - a*ny
      Rv(4,1) =  H - vn*a

      Rv(1,2) = zero
      Rv(2,2) = a*tx
      Rv(3,2) = a*ty
      Rv(4,2) = vt*a

      Rv(1,3) = one
      Rv(2,3) = vx
      Rv(3,3) = vy
      Rv(4,3) = half*(vx*vx+vy*vy)

      Rv(1,4) = one
      Rv(2,4) = vx + a*nx
      Rv(3,4) = vy + a*ny
      Rv(4,4) =  H + vn*a

      !Dissipation Term
      diss = zero
      do i=1,4
      do j=1,4
      diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
      end do
      end do

      !Compute the flux.
      fL(1) = rhoL*vnL
      fL(2) = rhoL*vnL * vxL + pL*nx
      fL(3) = rhoL*vnL * vyL + pL*ny
      fL(4) = rhoL*vnL *  HL

      fR(1) = rhoR*vnR
      fR(2) = rhoR*vnR * vxR + pR*nx
      fR(3) = rhoR*vnR * vyR + pR*ny
      fR(4) = rhoR*vnR *  HR

      roe = half * (fL + fR - diss)

 end function roe

!*****************************************************************************
!* -- Rotated-Roe-HLL Flux Function ---
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* The most robust Riemann solver known to the author (in terms of nonlinear
!* instability such as carbuncle).
!*
!* NB: At a boundary face, need to switch to a geometric normal vector:
!*               (nx2,ny2)=(nx, ny), (nx1,ny1)=(-ny,nx).
!*     This is not implemented here. It requires information on whether
!*     the geometric normal, (nx,ny), is on a boundary face or not.
!*     It shouldn't be difficult for you to implement it.
!*
!* Katate Masatsuka, February 2010. http://www.cfdbooks.com
!*****************************************************************************
 function rotated_rhll(uL, uR, nhat)
 double precision :: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
 double precision :: nhat(2)      !  Input: face normal vector
 double precision :: nx, ny          !  Input: face normal vector, [nx, ny] (Left-to-Right)
 double precision :: rotated_rhll(4) ! Output: Rotated_RHLL flux function.
!Local constants
 double precision :: gam                          ! Ratio of specific heat.
 double precision :: zero, fifth, half, one, two    ! Numbers
 double precision :: eps                            !
!Local variables
 double precision :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 double precision :: tx, ty                         ! Tangent vector (taken as n1)
 double precision :: alpha1, alpha2                 ! Projections of the new normals
 double precision :: vxL, vxR, vyL, vyR             ! Velocity components.
 double precision :: rhoL, rhoR, pL, pR             ! Primitive variables.
 double precision :: vnL, vnR, vtL, vtR             ! Normal and tagent velocities
 double precision :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
 double precision :: RT,rho,vx,vy,H,a               ! Roe-averages
 double precision :: vn, vt                         ! Normal and tagent velocities(Roe-average)
 double precision :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 double precision :: abs_dq                         ! Magnitude of the velocity difference
 double precision :: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigevectors
 double precision :: SRp,SLm                        ! Wave speeds for the HLL part
 double precision :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 double precision :: temp
 integer :: i, j

 ! Separating normal vector
 nx = nhat(1)
 ny = nhat(2)

!Constants.
     gam = 1.4
      zero = 0.0
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0
       eps = 1.0e-5 ! 1.0e-12 in the original paper (double precision)

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
      pL = (gam-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gam*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
      pR = (gam-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gam*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

     vnL = vxL*nx + vyL*ny
     vnR = vxR*nx + vyR*ny

!Compute the flux.
   fL(1) = rhoL*vnL
   fL(2) = rhoL*vnL * vxL + pL*nx
   fL(3) = rhoL*vnL * vyL + pL*ny
   fL(4) = rhoL*vnL *  HL

   fR(1) = rhoR*vnR
   fR(2) = rhoR*vnR * vxR + pR*nx
   fR(3) = rhoR*vnR * vyR + pR*ny
   fR(4) = rhoR*vnR *  HR

!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during
!     a steady calculation to fully make it converge. For time-accurate
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  if ( abs_dq > eps) then
       nx1 = (vxR-vxL)/abs_dq
       ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny
    ny1 =  nx
  endif
    alpha1 = nx * nx1 + ny * ny1
!   To make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! Take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
!   To make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
    alpha2 = temp * alpha2

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
     RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
     vx = (vxL+RT*vxR)/(one+RT)
     vy = (vyL+RT*vyR)/(one+RT)
      H = ( HL+RT* HR)/(one+RT)
      a = sqrt( (gam-one)*(H-half*(vx*vx+vy*vy)) )
     vn = vx*nx2+vy*ny2
     vt = vx*nx1+vy*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2
    vnR = vxR*nx2 + vyR*ny2
    vtL = vxL*nx1 + vyL*ny1
    vtR = vxR*nx1 + vyR*ny1

   drho = rhoR - rhoL
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
  abs_ws  = abs(ws)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, vtR + aR, vt + a)
   SLm = min( zero, vtL - aL, vt - a)

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1

  Rv(1,1) = one
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the Rotated-RHLL flux.
  rotated_rhll = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

 end function rotated_rhll

end module riemann
