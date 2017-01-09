module linalg
  use mesh_class, only : mesh
  implicit none
  private
  public :: thomas,norml2

  contains

    !------------------------------------------------------
    ! Subroutine for solving a tridiagonal matrix using
    ! Thomas algorithm.
    !
    ! IMPORTANT: You must allocate b to be the same size
    ! as d.  If you allocate b(2:n), the subroutine doesn't
    ! work properly.  I'm not sure what happens, but it
    ! breaks.
    !
    ! NOTES: Diagonal is d.  b is just below diagonal,
    !        a is above the diagonal, and c is the rhs.
    !
    !  [ d1 a1  0    ... 0] {x1}   {c1}
    !  [ b2 d2 a2 0  ... 0] {x2}   {c2}
    !  [  0 b3 d3 a3 ... 0] { .} = { .}
    !  [ .  .  .     ...am] { .}   { .}
    !  [ .  .  .     ...dn] { .}   { .}
    !------------------------------------------------------
    subroutine thomas(a,b,c,d,x)
      implicit none
      double precision, dimension(:), intent(in)    :: a,b
      double precision, dimension(:), intent(inout) :: c,d
      double precision, dimension(:), intent(out)   :: x
      integer :: j,nj

      ! Finding size of array
      nj = size(d)

      ! Converting tridiagonal system to upper triangular
      do j=2,nj
        d(j) = d(j) - b(j)/d(j-1)*a(j-1)
        c(j) = c(j) - b(j)/d(j-1)*c(j-1)
      end do

      ! Using back substitution to find the solution
      x(nj) = c(nj)/d(nj)
      do j=nj-1,1,-1
        x(j) = (c(j) - a(j)*x(j+1))/d(j)
      end do

    end subroutine thomas

    !------------------------------------------------------
    ! Function for computing the l2 norm of a field.
    ! This function is used in solvers.f90 to compute the
    ! norm of the residual only.  If it is going to be
    ! used for something else, it should be modified.
    !------------------------------------------------------
    pure function norml2(grid,r) result(n)
      implicit none
      type(mesh),       intent(in)                 :: grid
      double precision, intent(in), dimension(:,:) :: r
      double precision                             :: n
      integer                                      :: i,j

      ! Initializing norm to zero
      n = 0.0d0

      ! Computing norm
      do j=1,grid%nelemj
        do i=1,grid%nelemi
          n = n + r(i,j)**2*grid%elem(i,j)%area**3   ! r includes area already
        end do
      end do

      ! Computing the square root of the norm
      n = sqrt(n)

    end function norml2

end module linalg
