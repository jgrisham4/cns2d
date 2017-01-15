test_suite utils

! Global variables
double precision, parameter :: ptx = 0.528761d0
double precision, parameter :: pty = 0.818584d0
double precision            :: p4(2)

setup
end setup

teardown
end teardown

! Testing function for reflecting coordinate about line
test reflection
  double precision, dimension(2) :: p1,p2,p3

  ! Setting coordinates for points
  p1(1) = 0.0d0
  p1(2) = 1.0d0
  p2(1) = 1.1d0
  p2(2) = 1.45d0
  p3(1) = 0.25
  p3(2) = 1.5
  p4 = reflect(p1,p2,p3)

  ! Checking equality
  assert_equal_within(p4(1),ptx,1.0e-6)
  assert_equal_within(p4(2),pty,1.0e-6)
end test

! Testing function which finds the max states
test elem_max
  double precision, allocatable :: u(:,:,:),umax_comp(:,:,:)
  double precision              :: umax = 100.0d0

  ! Defining array
  allocate(u(3,2,4))
  allocate(umax_comp(3,2,4))
  u(1,1,:) = (/ -5.0d0, -0.5d0, 30.0d0, 0.0d0 /)
  u(2,1,:) = (/ -4.0d0, -2.0d0, 20.0d0, 0.0d0 /)
  u(3,1,:) = (/ -3.0d0, umax, 10.0d0, 0.0d0   /)
  u(1,2,:) = (/ -5.0d0, 1.0d0, 30.0d0, 0.0d0  /)
  u(2,2,:) = (/ -4.0d0, 1.0d0, 20.0d0, 0.0d0  /)
  u(3,2,:) = (/ -3.0d0, 1.0d0, 10.0d0, 0.0d0  /)

  ! Computing max states
  umax_comp = compute_elem_max(u,3,2)

  ! Checking equality
  assert_equal_within(umax_comp(3,1,2),umax,1.0e-10)

end test elem_max

! Testing function which finds the min states
test elem_min
  double precision, allocatable :: u(:,:,:),umin_comp(:,:,:)
  double precision              :: umin = -100.0d0

  ! Defining array
  allocate(u(3,2,4))
  allocate(umin_comp(3,2,4))
  u(1,1,:) = (/ -5.0d0, -0.5d0, 30.0d0, 0.0d0 /)
  u(2,1,:) = (/ -4.0d0, -2.0d0, 20.0d0, 0.0d0 /)
  u(3,1,:) = (/ -3.0d0, umin, 10.0d0, 0.0d0   /)
  u(1,2,:) = (/ -5.0d0, 1.0d0, 30.0d0, 0.0d0  /)
  u(2,2,:) = (/ -4.0d0, 1.0d0, 20.0d0, 0.0d0  /)
  u(3,2,:) = (/ -3.0d0, 1.0d0, 10.0d0, 0.0d0  /)

  ! Computing min states
  umin_comp = compute_elem_min(u,3,2)

  ! Checking equality
  assert_equal_within(umin_comp(3,1,2),umin,1.0e-10)

end test elem_min

end test_suite
