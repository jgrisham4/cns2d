program test_ghost
  use mesh_class, only : mesh, read_from_file, write_to_tec, preprocess
  implicit none
  character (len=30) :: mesh_file,tec_file
  integer            :: i,j

  ! Creating mesh object
  type(mesh) :: grid

  ! Defining inputs
  !mesh_file = "simple.cgns"
  !tec_file  = "simple.tec"
  mesh_file = "duct.cgns"
  tec_file  = "duct.tec"
  
  ! Reading mesh
  call read_from_file(grid,mesh_file)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Writing ghost cells to tecplot files
  open(2,file='ghost_bottom.tec')
  write(2,*) 'variables=x,y'
  write(2,'(a,i4,a)') 'zone i=',grid%imax+4,' j=2'
  do j=-1,0
    do i=-1,grid%imax+2
      write(2,'(es14.6,a,es14.6)') grid%x(i,j),' ',grid%y(i,j)
    end do
  end do
  close(2)

  open(3,file='ghost_top.tec')
  write(3,*) 'variables=x,y'
  write(3,'(a,i4,a)') 'zone i=',grid%imax+4,' j=2'
  do j=grid%jmax+1,grid%jmax+2
    do i=-1,grid%imax+2
      write(3,'(es14.6,a,es14.6)') grid%x(i,j),' ',grid%y(i,j)
    end do
  end do
  close(3)

  open(2,file='ghost_left.tec')
  write(2,*) 'variables=x,y'
  write(2,'(a,i4)') 'zone i=2 j=', grid%jmax
  do j=1,grid%jmax
    do i=-1,0
      write(2,'(es14.6,a,es14.6)') grid%x(i,j),' ',grid%y(i,j)
    end do
  end do
  close(2)

  open(3,file='ghost_right.tec')
  write(3,*) 'variables=x,y'
  write(3,'(a,i4)') 'zone i=2 j=', grid%jmax
  do j=1,grid%jmax
    do i=grid%imax+1,grid%imax+2
      write(3,'(es14.6,a,es14.6)') grid%x(i,j),' ',grid%y(i,j)
    end do
  end do
  close(3)

  ! Writing mesh to file
  call write_to_tec(grid,tec_file)

end program 
