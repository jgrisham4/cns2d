test_suite grad

! Global variables
double precision, parameter :: pi = 4.0d0*atan(1.0d0)

! Setting up (executed before each test)
setup
end setup

! Tearing down (executed after each test)
teardown
end teardown

! Testing gradient calculation on uniform grid
test gradient
  use mesh_class, only : mesh,read_from_file,preprocess
  implicit none

  type(mesh), dimension(3)         :: grids
  double precision, allocatable    :: f(:,:,:)
  double precision, allocatable    :: gradf(:,:,:)
  double precision                 :: nerrx(3),nerry(3),dx(3)
  double precision, dimension(3)   :: xv,yv
  double precision                 :: A(3,2),xstar(2),b(3)
  double precision                 :: Am(2,2),bm(2)
  integer                          :: i,j,k,aerr
  character (len=30), dimension(3) :: fnames
  character (len=30), dimension(3) :: tecnames

  ! Setting file names
  fnames(1)   = "../tests/square050.cgns"
  fnames(2)   = "../tests/square100.cgns"
  fnames(3)   = "../tests/square200.cgns"
  tecnames(1) = "../tests/square050grad.tec"
  tecnames(2) = "../tests/square100grad.tec"
  tecnames(3) = "../tests/square200grad.tec"

  ! Computing dx
  dx(1) = 1.0d0/50.0d0
  dx(2) = 1.0d0/100.0d0
  dx(3) = 1.0d0/200.0d0
  
  ! Looping over meshes
  do k=1,3

    ! Reading mesh
    call read_from_file(grids(k),fnames(k))

    ! Preprocessing mesh
    call preprocess(grids(k))

    ! Creating data to take gradient of
    if (allocated(f)) deallocate(f)
    if (allocated(gradf)) deallocate(gradf)
    allocate(f(grids(k)%nelemi,grids(k)%nelemj,4),stat=aerr)
    if (aerr.ne.0) then
      print *, "Error: Can't allocate memory for f."
      stop
    end if
    allocate(gradf(grids(k)%nelemi,grids(k)%nelemj,8),stat=aerr)
    if (aerr.ne.0) then
      print *, "Error: Can't allocate memory for gradf."
      stop
    end if
    do j=1,grids(k)%nelemj
      do i=1,grids(k)%nelemi
        f(i,j,1) = sin(2.0d0*pi*grids(k)%elem(i,j)%xc)*sin(2.0d0*pi*grids(k)%elem(i,j)%yc)
        f(i,j,2) = 0.0d0
        f(i,j,3) = 0.0d0
        f(i,j,4) = 0.0d0
      end do
    end do

    ! Computing gradient
    call compute_gradient(grids(k),f,gradf)

    ! Computing norm of error in x-direction
    nerrx(k) = 0.0d0
    nerry(k) = 0.0d0
    do j=1,grids(k)%nelemj
      do i=1,grids(k)%nelemi
        nerrx(k) = nerrx(k) + (gradf(i,j,1) - 2.0d0*pi*cos(2.0*pi*grids(k)%elem(i,j)%xc)*sin(2.0*pi*grids(k)%elem(i,j)%yc))**2*grids(k)%elem(i,j)%area
        nerry(k) = nerry(k) + (gradf(i,j,5) - 2.0d0*pi*sin(2.0*pi*grids(k)%elem(i,j)%xc)*cos(2.0*pi*grids(k)%elem(i,j)%yc))**2*grids(k)%elem(i,j)%area
      end do
    end do
    nerrx(k) = sqrt(nerrx(k))
    nerry(k) = sqrt(nerry(k))
    !print *, "Error in df/dx = ", nerrx(k)

    ! Writing result to file
    open(2,file=tecnames(k))
    write(2,'(a)') 'title="gradient data"'
    write(2,'(a)') 'variables="x","y","f","f_x","f_y"'
    write(2,'(a,i5,a,i5)') 'zone i=', grids(k)%imax, ' j=', grids(k)%jmax
    write(2,'(a)') 'datapacking=block'
    write(2,'(a)') 'varlocation=([3,4,5]=cellcentered)'
    do j=1,grids(k)%jmax
      do i=1,grids(k)%imax
        write(2,'(es25.10)',advance='no') grids(k)%x(i,j)
      end do
      write(2,'(a)') " "
    end do
    do j=1,grids(k)%jmax
      do i=1,grids(k)%imax
        write(2,'(es25.10)',advance='no') grids(k)%y(i,j)
      end do
      write(2,'(a)') " "
    end do
    do j=1,grids(k)%nelemj
      do i=1,grids(k)%nelemi
        write(2,'(es25.10)',advance='no') f(i,j,1)
      end do
      write(2,'(a)') " "
    end do
    do j=1,grids(k)%nelemj
      do i=1,grids(k)%nelemi
        write(2,'(es25.10)',advance='no') gradf(i,j,1)
      end do
      write(2,'(a)') " "
    end do
    do j=1,grids(k)%nelemj
      do i=1,grids(k)%nelemi
        write(2,'(es25.10)',advance='no') gradf(i,j,5)
      end do
      write(2,'(a)') " "
    end do

  end do

  ! Closing file
  close(2)

  ! Computing the order-of-accuracy for the x-component of the gradient
  xv = log10(dx)
  yv = log10(nerrx)
  A(1,1) = xv(1)
  A(1,2) = 1.0d0
  A(2,1) = xv(2)
  A(2,2) = 1.0d0
  A(3,1) = xv(3)
  A(3,2) = 1.0d0
  b = yv
  Am = matmul(transpose(A),A)
  bm = matmul(transpose(A),b)
  xstar(1) = (bm(1)*Am(2,2) - Am(1,2)*bm(2))/(Am(1,1)*Am(2,2) - Am(1,2)*Am(2,1))
  print *, "order of accuracy for x-component of gradient is ", xstar(1)
  assert_equal_within(xstar(1),2.0d0,1.0e-1)

  ! Computing the order-of-accuracy for the y-component of the gradient
  xv = log10(dx)
  yv = log10(nerry)
  b = yv
  bm = matmul(transpose(A),b)
  xstar(1) = (bm(1)*Am(2,2) - Am(1,2)*bm(2))/(Am(1,1)*Am(2,2) - Am(1,2)*Am(2,1))
  print *, "order of accuracy for y-component of gradient is ", xstar(1)
  assert_equal_within(xstar(1),2.0d0,1.0e-1)

end test

end test_suite
