program test_cgns

  use cgns
  double precision :: x(3,3), y(3,3)
  integer(kind=4) :: isize(2,3)
  character (len=30) :: basename,zonename
  integer :: imax,jmax

  ! Creating grid points for this simple example
  imax = 3
  jmax = 3
  do j=1,jmax
    do i=1,imax
      x(i,j) = dble(i-1)
      y(i,j) = dble(j-1)
    end do
  end do
  print *, 'Created simple 2-D grid points.'

  ! Opening a CGNS file
  call cg_open_f('simple.cgns',CG_MODE_WRITE,index_file,ier)
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
  call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateX',x,index_coord,ier)
  call cg_coord_write_f(index_file,index_base,index_zone,RealDouble,'CoordinateY',y,index_coord,ier)

  ! Closing cgns file
  call cg_close_f(index_file,ier)
  print *, 'Done writing mesh.'

end program test_cgns

