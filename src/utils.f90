!===========================================================
! This module contains some utility functions
!
! NOTE: The assumption of a calorically perfect gas is
!       made in the below conversions to and from
!       conservative and primitive variable vectors.
!===========================================================

module utils
  implicit none
  private
  public::w_to_u,u_to_w,reflect,nvec,compute_elem_max,compute_elem_min

  contains

    !------------------------------------------------------
    ! Function for converting a vector of primitive
    ! variables to a vector of conserved variables.
    ! Assuming a calorically perfect gas.
    !------------------------------------------------------
    function w_to_u(w,g) result(u)
      implicit none
      double precision, intent(in) :: w(4),g
      double precision             :: u(4)
      u(1) = w(1)
      u(2) = w(1)*w(2)
      u(3) = w(1)*w(3)
      u(4) = w(4)/(g-1.0d0)+0.5d0*w(1)*(w(2)**2 + w(3)**2)
    end function w_to_u

    !------------------------------------------------------
    ! Function for converting from vector of conserved
    ! variables to a vector of primitive variables.
    ! Assuming a calorically perfect gas.
    !------------------------------------------------------
    function u_to_w(u,g) result(w)
      implicit none
      double precision, intent(in) :: u(4),g
      double precision             :: w(4)
      w(1) = u(1)
      w(2) = u(2)/u(1)
      w(3) = u(3)/u(1)
      w(4) = (g-1.0d0)*(u(4)-0.5d0*w(1)*(w(2)**2 + w(3)**2))
    end function u_to_w

    !------------------------------------------------------
    ! Function for reflecting a vector about another vector
    ! pt1 and pt2 create a line about which pt3 is to be
    ! reflected
    !------------------------------------------------------
    pure function reflect(pt1,pt2,pt3) result(ptref)
      implicit none
      double precision, intent(in)  :: pt1(2)
      double precision, intent(in)  :: pt2(2)
      double precision, intent(in)  :: pt3(2)
      double precision              :: ptref(2)
      double precision,dimension(2) :: u,v,x

      ! Creating vectors
      u = pt3-pt1
      v = pt2-pt1

      ! Calculations
      x = u - (dot_product(u,v)/(norm2(v)**2))*v
      ptref = u - 2.0d0*x + pt1

    end function reflect

    !------------------------------------------------------
    ! Function for computing a normal vector given
    ! components of a vector
    !------------------------------------------------------
    pure function nvec(xc,yc) result(nhat)
      implicit none
      double precision, intent(in) :: xc      ! x-component
      double precision, intent(in) :: yc      ! y-component
      double precision             :: mag     ! magnitude of vector
      double precision             :: nhat(2)

      ! Computing magnitude of vector
      mag = sqrt(xc*xc + yc*yc)
      nhat(1) = xc/mag
      nhat(2) = yc/mag

    end function nvec

    !------------------------------------------------------
    ! Function for determining the max state for
    ! each element and the surrounding elements
    !------------------------------------------------------
    pure function compute_elem_max(u,ni,nj) result(um)
      implicit none
      double precision, allocatable, intent(in) :: u(:,:,:)
      double precision, allocatable             :: um(:,:,:)
      integer, intent(in)                       :: ni,nj
      integer                                   :: i,j,k,aer

      ! Allocating memory for um
      allocate(um(ni,nj,4),stat=aer)

      ! Looping over all elements
      do j=1,nj
        do i=1,ni
          if ((i.ne.1).and.(i.ne.ni).and.(j.ne.1).and.(j.ne.nj)) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if ((i.eq.1).and.(j.eq.1)) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i,j+1,k))
            end do
          else if ((i.eq.1).and.(j.eq.nj)) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i,j-1,k))
            end do
          else if ((i.eq.ni).and.(j.eq.1)) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i,j+1,k),u(i-1,j,k))
            end do
          else if ((i.eq.ni).and.(j.eq.nj)) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if (i.eq.1) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i,j-1,k))
            end do
          else if (i.eq.ni) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i,j+1,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if (j.eq.1) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i-1,j,k))
            end do
          else if (j.eq.nj) then
            do k=1,4
              um(i,j,k) = max(u(i,j,k),u(i+1,j,k),u(i-1,j,k),u(i,j-1,k))
            end do
          end if
        end do
      end do

    end function compute_elem_max

    !------------------------------------------------------
    ! Function for determining the min state for
    ! each element and the surrounding elements
    !------------------------------------------------------
    pure function compute_elem_min(u,ni,nj) result(um)
      implicit none
      double precision, allocatable, intent(in) :: u(:,:,:)
      double precision, allocatable             :: um(:,:,:)
      integer, intent(in)                       :: ni,nj
      integer                                   :: i,j,k,aer

      ! Allocating memory for um
      allocate(um(ni,nj,4),stat=aer)

      ! Looping over all elements
      do j=1,nj
        do i=1,ni
          if ((i.ne.1).and.(i.ne.ni).and.(j.ne.1).and.(j.ne.nj)) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if ((i.eq.1).and.(j.eq.1)) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i,j+1,k))
            end do
          else if ((i.eq.1).and.(j.eq.nj)) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i,j-1,k))
            end do
          else if ((i.eq.ni).and.(j.eq.1)) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i,j+1,k),u(i-1,j,k))
            end do
          else if ((i.eq.ni).and.(j.eq.nj)) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if (i.eq.1) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i,j-1,k))
            end do
          else if (i.eq.ni) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i,j+1,k),u(i-1,j,k),u(i,j-1,k))
            end do
          else if (j.eq.1) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i,j+1,k),u(i-1,j,k))
            end do
          else if (j.eq.nj) then
            do k=1,4
              um(i,j,k) = min(u(i,j,k),u(i+1,j,k),u(i-1,j,k),u(i,j-1,k))
            end do
          end if
        end do
      end do

    end function compute_elem_min

end module utils
