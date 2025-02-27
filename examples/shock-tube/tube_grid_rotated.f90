program test_cgns
  use cgns
  double precision, parameter   :: pi = 4.0d0*atan(1.0d0)
  double precision, allocatable :: x(:,:),y(:,:),xp(:,:),yp(:,:)
  double precision              :: xmin,xmax,ymin,ymax
  double precision              :: dx,dy
  double precision              :: theta
  integer(kind=4)               :: isize(2,3)
  character (len=30)            :: basename,zonename,fname
  integer :: imax,jmax,i,j,icelldim,iphysdim,index_base,index_file,index_zone,index_coord,ier,aerr

  ! Setting grid inputs
  imax = 201
  jmax = 21
  xmin = 0.0d0
  xmax = 1.0d0
  ymin = 0.0d0
  ymax = 0.2d0
  theta = 45.0d0*pi/180.0d0
  fname = "tube_rotated.cgns"

  ! Allocating memory for arrays
  allocate(x(imax,jmax),stat=aerr)
  if (aerr.ne.0) then
    print *, "Error: Can't allocate memory for x."
    stop
  end if
  allocate(y(imax,jmax),stat=aerr)
  if (aerr.ne.0) then
    print *, "Error: Can't allocate memory for y."
    stop
  end if
  allocate(xp(imax,jmax))
  allocate(yp(imax,jmax))

  ! Creating grid points for this simple example
  dx = (xmax - xmin)/dble(imax-1)
  dy = (ymax - ymin)/dble(jmax-1)
  do j=1,jmax
    do i=1,imax
      x(i,j) = dx*(dble(i-1)) + xmin
      y(i,j) = dy*(dble(j-1)) + ymin
    end do
  end do
  print *, 'Done creating 2-D grid points.'

  ! Transforming using the rotation matrix
  do j=1,jmax
    do i=1,imax
      xp(i,j) = x(i,j)*cos(theta) - y(i,j)*sin(theta)
      yp(i,j) = x(i,j)*sin(theta) + y(i,j)*cos(theta)
    end do
  end do

  ! Opening a CGNS file
  call cg_open_f(fname,CG_MODE_WRITE,index_file,ier)
  if (ier.ne.CG_OK) call cg_error_exit_f

  ! Creating base
  basename='Base'
  icelldim=2
  iphysdim=2
  call cg_base_write_f(index_file,basename,icelldim,iphysdim,index_base,ier)

  ! Defining a zone
  zonename = 'Zone 1'

  ! Setting the number of nodes in each direction
  isize(1,1) = imax
  isize(2,1) = jmax

  ! Setting the number of cells in each direction
  isize(1,2) = imax-1
  isize(2,2) = jmax-1

  ! Setting boundary vertex size (0 for structured)
  isize(1,3) = 0
  isize(2,3) = 0

  ! Creating a zone
  call cg_zone_write_f(index_file,index_base,zonename,isize,Structured,index_zone,ier)
  if (ier.ne.CG_OK) call cg_error_exit_f

  ! Writing data to zone
  call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateX',xp,index_coord,ier)
  call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateY',yp,index_coord,ier)

  ! Closing cgns file
  call cg_close_f(index_file,ier)
  print *, 'Done writing mesh.'

end program test_cgns

