test_suite linalg

! Setup portion
setup
end setup

! Teardown
teardown
end teardown

!-------------------------------------------------------------------------------
! Testing the thomas algorithm
!-------------------------------------------------------------------------------
test thomas_algorithm
  double precision, dimension(4) :: ud
  double precision               :: ld(5)
  double precision, dimension(5) :: d,rhs,sol,sol_act

  ! Filling in values on diagonals
  d(1) = 5.0d0
  d(2) = 4.0d0
  d(3) = 2.0d0
  d(4) = 1.0d0
  d(5) = 8.0d0
  ld(2) = 1.0d0
  ld(3) = 5.0d0
  ld(4) = -10.0d0
  ld(5) = 2.0d0
  ud(1) = -2.0d0
  ud(2) = 5.0d0
  ud(3) = -10.0d0
  ud(4) = 1.0d0

  ! Setting the actual solution
  sol_act(:) = 1.0d0

  ! Setting RHS
  rhs(1) = 3.0d0
  rhs(2) = 10.0d0
  rhs(3) = -3.0d0
  rhs(4) = -8.0d0
  rhs(5) = 10.0d0

  ! Computing the solution using thomas algorithm
  call thomas(ud,ld,rhs,d,sol)

  ! Making sure the solution is correct
  assert_equal_within(sol(1),sol_act(1),1.0e-10)
  assert_equal_within(sol(2),sol_act(2),1.0e-10)
  assert_equal_within(sol(3),sol_act(3),1.0e-10)
  assert_equal_within(sol(4),sol_act(4),1.0e-10)
  assert_equal_within(sol(5),sol_act(5),1.0e-10)

end test

!-------------------------------------------------------------------------------
! Testing the residual calculation
! For this case, I am taking the norm of ones over the entire domain.  This
! should return the area of the domain.
!-------------------------------------------------------------------------------
test residual_norm
  use mesh_class, only : mesh,read_from_file,preprocess
  character (len=30) :: grid_file
  double precision, allocatable :: resid(:,:)
  double precision :: n
  type(mesh) :: grid
  integer :: i,j

  ! Reading grid from file
  grid_file = "../tests/square050.cgns"
  call read_from_file(grid,grid_file)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Allocating memory for fake residual
  allocate(resid(grid%nelemi,grid%nelemj))

  ! Defining ones at the cell center for each element
  do j=1,grid%nelemj
    do i=1,grid%nelemi
      resid(i,j) = 1.0d0/grid%elem(i,j)%area
    end do
  end do

  ! Computing the norm
  n = norml2(grid,resid)

  ! Making sure that the result is approximately equal to the area of the domain
  assert_equal_within(n,1.0d0,1.0e-10)

end test

end test_suite
