program test_reflection
  use utils, only : reflect
  implicit none
  double precision, dimension(2) :: p1,p2,p3,p4

  p1(1) = 0.0d0
  p1(2) = 1.0d0
  p2(1) = 1.0d0
  p2(2) = 1.25d0
  p3(1) = 0.25d0
  p3(2) = 1.5d0
  p4 = reflect(p1,p2,p3)
  
  print *, 'p4 = ', p4

end program
