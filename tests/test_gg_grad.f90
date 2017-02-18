program test_gg_grad
  use mesh_class, only : mesh,read_from_file,preprocess,edge
  use utils,      only : w_to_u
  use grad,       only : compute_face_gradients
  implicit none

  type(mesh)         :: grid
  type(edge)         :: edgetmp
  integer            :: i,j,k,aerr
  character (len=30) :: fname
  character (len=30) :: v_tecname,h_tecname
  double precision   :: xtmp,ytmp

  ! Defining pi
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)

  ! Parameters for MMS
  double precision, parameter :: L    = 1.0d0
  double precision, parameter :: r0   = 1.0d0
  double precision, parameter :: u0   = 800.0d0
  double precision, parameter :: v0   = 800.0d0
  double precision, parameter :: p0   = 100000.0d0
  double precision, parameter :: rx   = 0.15d0
  double precision, parameter :: ux   = 50.0d0
  double precision, parameter :: vx   = -75.0d0
  double precision, parameter :: px   = 2000.0d0
  double precision, parameter :: ry   = -0.1d0
  double precision, parameter :: uy   = -30.0d0
  double precision, parameter :: vy   = 40.0d0
  double precision, parameter :: py   = 5000.0d0
  double precision, parameter :: arx  = 1.0d0
  double precision, parameter :: aux  = 1.5d0
  double precision, parameter :: avx  = 0.5d0
  double precision, parameter :: apx  = 2.0d0/3.0d0
  double precision, parameter :: ary  = 0.5d0
  double precision, parameter :: auy  = 0.6d0
  double precision, parameter :: avy  = 1.5d0
  double precision, parameter :: apy  = 1.0d0

  ! Gas properties
  double precision, parameter :: g    = 1.4d0
  double precision, parameter :: Rgas = 287.0d0

  ! Reading mesh
  fname = "../tests/square100.cgns"
  call read_from_file(grid,fname)

  ! Preprocessing mesh
  call preprocess(grid)

  ! Setting data using manufactured solution
  do j=0,grid%nelemj+1
    do i=0,grid%nelemi+1

      ! Getting coordinates of the cell center
      xtmp = grid%elem(i,j)%xc
      ytmp = grid%elem(i,j)%yc

      ! Setting primitive variables at the cell center
      grid%elem(i,j)%w(1) = r0 + rx*sin(arx*pi*xtmp/L) + ry*cos(ary*pi*ytmp/L)
      grid%elem(i,j)%w(2) = u0 + ux*sin(aux*pi*xtmp/L) + uy*cos(auy*pi*ytmp/L)
      grid%elem(i,j)%w(3) = v0 + vx*sin(avx*pi*xtmp/L) + vy*cos(avy*pi*ytmp/L)
      grid%elem(i,j)%w(4) = p0 + px*sin(apx*pi*xtmp/L) + py*cos(apy*pi*ytmp/L)

      ! Converting to conserved variables
      grid%elem(i,j)%u = w_to_u(grid%elem(i,j)%w,g)

    end do
  end do

  ! Computing the gradient
  call compute_face_gradients(grid,Rgas,g)

  ! Writing the results to tecplot file
  v_tecname = "grad_v.tec"
  h_tecname = "grad_h.tec"
  open(2,file=v_tecname)
  write(2,'(a)') 'title="vertical face gradients"'
  write(2,'(a)') 'variables="x","y","dudx","dudy","dvdx","dvdy","dudx_e","dudy_e"'
  write(2,'(a,i5,a,i5,a)') 'zone i=', grid%imax, ' j=', grid%nelemj, ' f=point'
  do j=1,grid%nelemj
    do i=1,grid%imax
      edgetmp = grid%edges_v(i,j)
      xtmp = edgetmp%xm
      ytmp = edgetmp%ym
      write(2,'(6(es25.10,1x))') edgetmp%xm, edgetmp%ym, edgetmp%gradu(1), edgetmp%gradu(2), &
        edgetmp%gradv(1), edgetmp%gradv(2), aux*pi*ux*cos(aux*pi*xtmp/L)/L, &
        -auy*pi*uy*sin(auy*pi*ytmp/L)/L
    end do
  end do

  close(2)

end program test_gg_grad
