test_suite utils

! Global variables
double precision, parameter :: ptx = 0.528761d0
double precision, parameter :: pty = 0.818584d0
double precision            :: p4(2)

setup
  double precision, dimension(2) :: p1,p2,p3

  ! Setting coordinates for points
  p1(1) = 0.0d0
  p1(2) = 1.0d0
  p2(1) = 1.1d0
  p2(2) = 1.45d0
  p3(1) = 0.25
  p3(2) = 1.5
  p4 = reflect(p1,p2,p3)

end setup

teardown
end teardown

test reflection
  assert_equal_within(p4(1),ptx,1.0e-6)
  assert_equal_within(p4(2),pty,1.0e-6)
end test

end test_suite
