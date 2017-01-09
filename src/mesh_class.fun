test_suite mesh_class

! Creating mesh object (global scope)
character (len=30) :: mesh_file
double precision, parameter :: pi = 3.1415927
double precision, parameter :: theta = 12.0d0

setup

  mesh_file = "../tests/simple.cgns"

end setup

teardown
end teardown

!---------------------------------------------------------------
! Test for simple CGNS mesh read
!---------------------------------------------------------------
test cgns_read
  type(mesh) :: grid

  ! Reading simple 3x3 mesh
  call read_from_file(grid,mesh_file)

  ! Checking coordinates
  assert_equal_within(grid%x(1,1),0.0d0,1.0e-10)
  assert_equal_within(grid%x(2,1),1.0d0,1.0e-10)
  assert_equal_within(grid%x(3,1),2.0d0,1.0e-10)
  assert_equal_within(grid%y(1,1),0.0d0,1.0e-10)
  assert_equal_within(grid%y(1,2),1.0d0,1.0e-10)
  assert_equal_within(grid%y(1,3),2.0d0,1.0e-10)

  ! Checking number of nodes
  assert_equal(grid%imax*grid%jmax,16)

  ! Checking number of elements
  assert_equal(grid%nelemi,3)
  assert_equal(grid%nelemj,3)

end test

!---------------------------------------------------------------
! Test for mesh preprocessing
!---------------------------------------------------------------
test mesh_preprocess
  type(mesh) :: grid

  ! Reading simple 3x3 mesh
  call read_from_file(grid,mesh_file)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Checking centroids
  assert_equal_within(grid%elem(1,1)%xc,0.5d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%yc,0.5d0,1.0e-10)
  assert_equal_within(grid%elem(2,1)%xc,1.5d0,1.0e-10)
  assert_equal_within(grid%elem(2,1)%yc,0.5d0,1.0e-10)

  ! Checking area
  assert_equal_within(grid%elem(1,1)%area,1.0d0,1.0e-10)

  ! Checking face normals
  assert_equal_within(grid%elem(1,1)%n(1,1),0.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(2,1),-1.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(1,2),1.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(2,2),0.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(1,3),0.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(2,3),1.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(1,4),-1.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,1)%n(2,4),0.0d0,1.0e-10)

  ! Checking metrics
  assert_equal_within(grid%elem(3,3)%dxdxi,1.0d0,1.0e-10)
  assert_equal_within(grid%elem(3,3)%dydxi,0.0d0,1.0e-10)
  assert_equal_within(grid%elem(3,3)%dxdeta,0.0d0,1.0e-10)
  assert_equal_within(grid%elem(3,3)%dydeta,1.0d0,1.0e-10)
  assert_equal_within(grid%elem(3,3)%detJ,1.0d0,1.0e-10)

end test

!---------------------------------------------------------------
! Test for mesh preprocessing on a rotated mesh
!---------------------------------------------------------------
test mesh_preprocess_rotated
  type(mesh) :: grid
  double precision :: th

  ! Computing theta in radians
  th = theta*pi/180.0d0

  ! Reading simple 3x3 mesh
  call read_from_file(grid,"../tests/rotated.cgns")

  ! Preprocessing mesh
  call preprocess(grid)

  ! Checking area
  assert_equal_within(grid%elem(1,1)%area,1.0d0,1.0e-10)
  assert_equal_within(grid%elem(1,grid%nelemj)%area,2.0d0,1.0e-10)

  ! Checking face normals
  assert_equal_within(grid%elem(1,1)%n(1,1), sin(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(2,1),-cos(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(1,2), cos(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(2,2), sin(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(1,3),-sin(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(2,3), cos(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(1,4),-cos(th),1.0e-8)
  assert_equal_within(grid%elem(1,1)%n(2,4),-sin(th),1.0e-8)

  ! Checking metrics
  assert_equal_within(grid%elem(1,grid%nelemj)%detJ,1.75d0,1.0e-10)

  ! Checking edge length
  assert_equal_within(grid%edges_v(1,grid%nelemj)%length,2.0d0,1.0e-10)

end test

end test_suite
