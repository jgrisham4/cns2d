!===========================================================
! This module contains the definition for the mesh class.
!===========================================================

module mesh_class
  use cgns
  use utils, only : reflect, nvec
  implicit none
  private
  public :: mesh,element,read_from_file,write_to_tec,preprocess

  !---------------------------------------------------------
  ! Edge class definition
  !---------------------------------------------------------
  type edge
    double precision :: length  ! length of interface
    double precision :: xm,ym   ! x- and y-coordinates of midpoint
    double precision :: uL(4)   ! Left interface state
    double precision :: uR(4)   ! Right interface state
    double precision :: flux(4) ! value of the flux at the interface (i.e., f.n)
  end type edge

  !---------------------------------------------------------
  ! Element class definition
  !---------------------------------------------------------
  type element
    double precision :: xc,yc    ! Coordinates of cell centroid
    double precision :: area     ! Area of the cell
    double precision :: n(2,4)   ! Face normals
    double precision :: dxdxi    ! dx/dxi
    double precision :: dydxi    ! dy/dxi
    double precision :: dxdeta   ! dx/deta
    double precision :: dydeta   ! dy/deta
    double precision :: detJ     ! Jacobian determinant
    double precision :: u(4)     ! Vector of conserved variables
    double precision :: u0(4)    ! Vector of conserved variables at previous timestep
    double precision :: w(4)     ! Vector of primitive variables
    double precision :: dudx(4)  ! Gradient of u in x-direction
    double precision :: dudy(4)  ! Gradient of u in y-direction
    double precision :: dTdx     ! Gradient of T in x-direction
    double precision :: dVdx(2)  ! Gradient of velocity in x-direction
    double precision :: dVdy(2)  ! Gradient of velocity in y-direction
    double precision :: dTdy     ! Gradient of T in y-direction
    double precision :: wi(4,4)  ! Primitive states at interfaces (1-bottom, 2-right,...)
  end type element

  !---------------------------------------------------------
  ! Mesh class definition
  !---------------------------------------------------------
  type mesh
    integer                       :: imax,jmax      ! Number of nodes in i and j
    integer                       :: nelemi,nelemj  ! Number of elements in i and j
    double precision, allocatable :: x(:,:),y(:,:)  ! Node coordinates
    type(element), allocatable    :: elem(:,:)      ! Array of element objects
    type(edge), allocatable       :: edges_h(:,:)   ! Horizontal edges
    type(edge), allocatable       :: edges_v(:,:)   ! Vertical edges
  end type mesh

  contains

    !-------------------------------------------------------
    ! Subroutine for reading the mesh from a CGNS file
    ! This only reads x,y values from CGNS file
    !-------------------------------------------------------
    subroutine read_from_file(this,file_name)
      implicit none
      type(mesh), intent(inout)     :: this
      character (len=*)             :: file_name
      character (len=30)            :: zonename
      integer                       :: aerr,i,j
      integer                       :: index_file,index_base,index_zone
      integer                       :: ier
      !integer(kind=8)               :: isize(2,3)          ! gfortran
      integer(kind=4)               :: isize(3,3)         ! ifort
      integer(kind=4)               :: irmin(3),irmax(3)
      double precision, allocatable :: xtmp(:,:),ytmp(:,:)

      print *, "Reading mesh from ", file_name, "..."

      ! Opening CGNS file
      call cg_open_f(file_name,CG_MODE_READ,index_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      index_base = 1
      index_zone = 1

      ! Getting zone size
      call cg_zone_read_f(index_file,index_base,index_zone,zonename,isize,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      ! Lower range index
      irmin(1) = 1
      irmin(2) = 1

      ! Upper range index
      irmax(1)  = isize(1,1)
      irmax(2)  = isize(2,1)
      this%imax = irmax(1)
      this%jmax = irmax(2)
      this%nelemi = this%imax-1
      this%nelemj = this%jmax-1

      ! Some debugging info
      !write (*,'(3i4)') irmin
      !write (*,'(3i4)') irmax
      !write (*,'(a)') "isize = "
      !do j=1,3
      !  do i=1,3
      !    write (*,'(i4)',advance='no') isize(i,j)
      !  end do
      !  write (*,'(a)') " "
      !end do

      ! Allocating memory for x and y
      allocate(xtmp(this%imax,this%jmax),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for xtmp."
        stop
      end if
      allocate(ytmp(this%imax,this%jmax),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for ytmp."
        stop
      end if
      allocate(this%x(-1:this%imax+2,-1:this%jmax+2),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for x."
        stop
      end if
      allocate(this%y(-1:this%imax+2,-1:this%jmax+2),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for y."
        stop
      end if

      ! Reading grid coordinates
      call cg_coord_read_f(index_file,index_base,index_zone,'CoordinateX',RealDouble,irmin,irmax,xtmp,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      call cg_coord_read_f(index_file,index_base,index_zone,'CoordinateY',RealDouble,irmin,irmax,ytmp,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f

      ! Closing CGNS file
      call cg_close_f(index_file,ier)
      if (ier.ne.CG_OK) call cg_error_exit_f
      print *, "Closed CGNS file."

      ! Copying grid points from xtmp and ytmp to x and y
      do j=1,this%jmax
        do i=1,this%imax
          this%x(i,j) = xtmp(i,j)
          this%y(i,j) = ytmp(i,j)
        end do
      end do

      ! Printing some info
      print *, "Done reading mesh."

    end subroutine read_from_file

    !-------------------------------------------------------
    ! Subroutine for writing the mesh to a Tecplot file
    !-------------------------------------------------------
    subroutine write_to_tec(this,file_name)
      implicit none
      type(mesh), intent(in) :: this
      character (len=*) :: file_name

      ! Declaring local variables
      integer :: i,j

      ! Opening file
      open(2,file=file_name)
      write(2,'(a)') 'variables=x,y'
      write(2,'(a,i5,a,i5)') 'zone i=',this%imax,' j=',this%jmax
      do j=1,this%jmax
        do i=1,this%imax
          write(2,'(es25.10,a,es25.10)') this%x(i,j),' ',this%y(i,j)
        end do
      end do
      close(2)

    end subroutine write_to_tec

    !-------------------------------------------------------
    ! Subroutine for preprocessing the mesh
    ! (i.e., computing geometric quantities and metrics)
    !-------------------------------------------------------
    subroutine preprocess(this)
      implicit none
      type(mesh), intent(inout)     :: this        ! mesh object
      integer                       :: aerr        ! allocate error indicator
      integer                       :: i,j         ! dummy indices
      double precision,dimension(2) :: p1,p2,p3,p4 ! Points used in reflections
      double precision              :: x1,x2,x3,x4 ! Coordinates of nodes for an element
      double precision              :: y1,y2,y3,y4 ! Coordinates of nodes for an element

      print *, "Preprocessing mesh..."

      ! Allocating memory
      allocate(this%elem(-1:(this%imax-1)+2,-1:(this%jmax-1)+2),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for elem."
        stop
      end if
      allocate(this%edges_h(this%nelemi,this%jmax),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for edges_h."
        stop
      end if
      allocate(this%edges_v(this%imax,this%nelemj),stat=aerr)
      if (aerr.ne.0) then
        print *, "Error: can't allocate memory for edges_v."
        stop
      end if

      ! Creating ghost nodes along bottom and top
      ! PARALLELIZE
      do i=1,this%imax-1

        ! Bottom first point
        p1(1) = this%x(i,1)
        p1(2) = this%y(i,1)
        p2(1) = this%x(i+1,1)
        p2(2) = this%y(i+1,1)
        p3(1) = this%x(i,2)
        p3(2) = this%y(i,2)
        p4 = reflect(p1,p2,p3)
        this%x(i,0) = p4(1)
        this%y(i,0) = p4(2)

        ! Bottom second point
        p3(1) = this%x(i,3)
        p3(2) = this%y(i,3)
        p4 = reflect(p1,p2,p3)
        this%x(i,-1) = p4(1)
        this%y(i,-1) = p4(2)

        ! Top first point
        p1(1) = this%x(i,this%jmax)
        p1(2) = this%y(i,this%jmax)
        p2(1) = this%x(i+1,this%jmax)
        p2(2) = this%y(i+1,this%jmax)
        p3(1) = this%x(i,this%jmax-1)
        p3(2) = this%y(i,this%jmax-1)
        p4 = reflect(p1,p2,p3)
        this%x(i,this%jmax+1) = p4(1)
        this%y(i,this%jmax+1) = p4(2)

        ! Top second point
        p3(1) = this%x(i,this%jmax-2)
        p3(2) = this%y(i,this%jmax-2)
        p4 = reflect(p1,p2,p3)
        this%x(i,this%jmax+2) = p4(1)
        this%y(i,this%jmax+2) = p4(2)

      end do

      ! Last two points on the bottom
      p1(1) = this%x(this%imax-1,1)
      p1(2) = this%y(this%imax-1,1)
      p2(1) = this%x(this%imax,1)
      p2(2) = this%y(this%imax,1)
      p3(1) = this%x(this%imax,2)
      p3(2) = this%y(this%imax,2)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax,0) = p4(1)
      this%y(this%imax,0) = p4(2)
      p3(1) = this%x(this%imax,3)
      p3(2) = this%y(this%imax,3)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax,-1) = p4(1)
      this%y(this%imax,-1) = p4(2)

      ! Last two points on the top
      p1(1) = this%x(this%imax-1,this%jmax)
      p1(2) = this%y(this%imax-1,this%jmax)
      p2(1) = this%x(this%imax,this%jmax)
      p2(2) = this%y(this%imax,this%jmax)
      p3(1) = this%x(this%imax,this%jmax-1)
      p3(2) = this%y(this%imax,this%jmax-1)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax,this%jmax+1) = p4(1)
      this%y(this%imax,this%jmax+1) = p4(2)
      p3(1) = this%x(this%imax,this%jmax-2)
      p3(2) = this%y(this%imax,this%jmax-2)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax,this%jmax+2) = p4(1)
      this%y(this%imax,this%jmax+2) = p4(2)

      ! Creating ghost nodes along left and right
      ! PARALLELIZE
      do j=1,this%jmax-1

        ! Left side first point
        p1(1) = this%x(1,j)
        p1(2) = this%y(1,j)
        p2(1) = this%x(1,j+1)
        p2(2) = this%y(1,j+1)
        p3(1) = this%x(2,j)
        p3(2) = this%y(2,j)
        p4 = reflect(p1,p2,p3)
        this%x(0,j) = p4(1)
        this%y(0,j) = p4(2)

        ! Left side second point
        p3(1) = this%x(3,j)
        p3(2) = this%y(3,j)
        p4 = reflect(p1,p2,p3)
        this%x(-1,j) = p4(1)
        this%y(-1,j) = p4(2)

        ! Right side first point
        p1(1) = this%x(this%imax,j)
        p1(2) = this%y(this%imax,j)
        p2(1) = this%x(this%imax,j+1)
        p2(2) = this%y(this%imax,j+1)
        p3(1) = this%x(this%imax-1,j)
        p3(2) = this%y(this%imax-1,j)
        p4 = reflect(p1,p2,p3)
        this%x(this%imax+1,j) = p4(1)
        this%y(this%imax+1,j) = p4(2)

        ! Right side second point
        p3(1) = this%x(this%imax-2,j)
        p3(2) = this%y(this%imax-2,j)
        p4 = reflect(p1,p2,p3)
        this%x(this%imax+2,j) = p4(1)
        this%y(this%imax+2,j) = p4(2)

      end do

      ! Last two points on the left
      p1(1) = this%x(1,this%jmax)
      p1(2) = this%y(1,this%jmax)
      p2(1) = this%x(1,this%jmax-1)
      p2(2) = this%y(1,this%jmax-1)
      p3(1) = this%x(2,this%jmax)
      p3(2) = this%y(2,this%jmax)
      p4 = reflect(p1,p2,p3)
      this%x(0,this%jmax) = p4(1)
      this%y(0,this%jmax) = p4(2)
      p3(1) = this%x(3,this%jmax)
      p3(2) = this%y(3,this%jmax)
      p4 = reflect(p1,p2,p3)
      this%x(-1,this%jmax) = p4(1)
      this%y(-1,this%jmax) = p4(2)

      ! Last two points on the right
      p1(1) = this%x(this%imax,this%jmax)
      p1(2) = this%y(this%imax,this%jmax)
      p2(1) = this%x(this%imax,this%jmax-1)
      p2(2) = this%y(this%imax,this%jmax-1)
      p3(1) = this%x(this%imax-1,this%jmax)
      p3(2) = this%y(this%imax-1,this%jmax)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax+1,this%jmax) = p4(1)
      this%y(this%imax+1,this%jmax) = p4(2)
      p3(1) = this%x(this%imax-2,this%jmax)
      p3(2) = this%y(this%imax-2,this%jmax)
      p4 = reflect(p1,p2,p3)
      this%x(this%imax+2,this%jmax) = p4(1)
      this%y(this%imax+2,this%jmax) = p4(2)

      ! Copying data to cells in corners
      ! Bottom left corner
      this%x(0,0)   = this%x(0,1)
      this%x(0,-1)  = this%x(0,1)
      this%x(-1,0)  = this%x(-1,1)
      this%x(-1,-1) = this%x(-1,1)
      this%y(0,0)   = this%y(1,0)
      this%y(-1,0)  = this%y(1,0)
      this%y(0,-1)  = this%y(1,-1)
      this%y(-1,-1) = this%y(1,-1)

      ! Bottom right corner
      this%x(this%imax+1,0)  = this%x(this%imax+1,1)
      this%x(this%imax+1,-1) = this%x(this%imax+1,1)
      this%x(this%imax+2,0)  = this%x(this%imax+2,1)
      this%x(this%imax+2,-1) = this%x(this%imax+2,1)
      this%y(this%imax+1,0)  = this%y(this%imax,0)
      this%y(this%imax+2,0)  = this%y(this%imax,0)
      this%y(this%imax+1,-1) = this%y(this%imax,-1)
      this%y(this%imax+2,-1) = this%y(this%imax,-1)

      ! Top right corner
      this%x(this%imax+1,this%jmax+1) = this%x(this%imax+1,this%jmax)
      this%x(this%imax+1,this%jmax+2) = this%x(this%imax+1,this%jmax)
      this%x(this%imax+2,this%jmax+1) = this%x(this%imax+2,this%jmax)
      this%x(this%imax+2,this%jmax+2) = this%x(this%imax+2,this%jmax)
      this%y(this%imax+1,this%jmax+1) = this%y(this%imax,this%jmax+1)
      this%y(this%imax+2,this%jmax+1) = this%y(this%imax,this%jmax+1)
      this%y(this%imax+1,this%jmax+2) = this%y(this%imax,this%jmax+2)
      this%y(this%imax+2,this%jmax+2) = this%y(this%imax,this%jmax+2)

      ! Top left corner
      this%x(0,this%jmax+1)  = this%x(0,this%jmax)
      this%x(0,this%jmax+2)  = this%x(0,this%jmax)
      this%x(-1,this%jmax+1) = this%x(-1,this%jmax)
      this%x(-1,this%jmax+2) = this%x(-1,this%jmax)
      this%y(0,this%jmax+1)  = this%y(1,this%jmax+1)
      this%y(-1,this%jmax+1) = this%y(1,this%jmax+1)
      this%y(0,this%jmax+2)  = this%y(1,this%jmax+2)
      this%y(-1,this%jmax+2) = this%y(1,this%jmax+2)

      ! Computing geometric quantities PARALLELIZE
      do j=-1,this%nelemj+2
        do i=-1,this%nelemi+2

          ! Getting coordinates of the nodes that belong to this element
          x1 = this%x(i,j)
          x2 = this%x(i+1,j)
          x3 = this%x(i+1,j+1)
          x4 = this%x(i,j+1)
          y1 = this%y(i,j)
          y2 = this%y(i+1,j)
          y3 = this%y(i+1,j+1)
          y4 = this%y(i,j+1)

          ! Constructing normal vectors for faces
          this%elem(i,j)%n(:,1) = nvec(y2-y1,x1-x2)
          this%elem(i,j)%n(:,2) = nvec(y3-y2,x2-x3)
          this%elem(i,j)%n(:,3) = nvec(y4-y3,x3-x4)
          this%elem(i,j)%n(:,4) = nvec(y1-y4,x4-x1)

          ! Computing centroids of elements
          this%elem(i,j)%xc = 0.25d0*(x1+x2+x3+x4)
          this%elem(i,j)%yc = 0.25d0*(y1+y2+y3+y4)

          ! Computing cell area
          this%elem(i,j)%area = 0.5d0*((x1-x3)*(y2-y4)+(x4-x2)*(y1-y3))

        end do
      end do

      ! Computing metrics on interior cells
      ! PARALLELIZE
      do j=2,this%nelemj-1
        do i=2,this%nelemi-1
          this%elem(i,j)%dxdxi  = 0.5d0*(this%elem(i+1,j)%xc - this%elem(i-1,j)%xc)
          this%elem(i,j)%dydxi  = 0.5d0*(this%elem(i+1,j)%yc - this%elem(i-1,j)%yc)
          this%elem(i,j)%dxdeta = 0.5d0*(this%elem(i,j+1)%xc - this%elem(i,j-1)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(this%elem(i,j+1)%yc - this%elem(i,j-1)%yc)
        end do
      end do

      ! Computing metrics on the vertical boundaries
      do j=1,this%nelemj

        ! left side
        i = 1
        this%elem(i,j)%dxdxi = 0.5d0*(-3.0d0*this%elem(i,j)%xc + 4.0d0*this%elem(i+1,j)%xc - this%elem(i+2,j)%xc)
        this%elem(i,j)%dydxi = 0.5d0*(-3.0d0*this%elem(i,j)%yc + 4.0d0*this%elem(i+1,j)%yc - this%elem(i+2,j)%yc)

        ! right size
        i = this%nelemi
        this%elem(i,j)%dxdxi = 0.5d0*(3.0d0*this%elem(i,j)%xc - 4.0d0*this%elem(i-1,j)%xc + this%elem(i-2,j)%xc)
        this%elem(i,j)%dydxi = 0.5d0*(3.0d0*this%elem(i,j)%yc - 4.0d0*this%elem(i-1,j)%yc + this%elem(i-2,j)%yc)

        ! Interior on vertical boundaries
        if (j.ne.1.and.j.ne.this%nelemj) then

          ! left side
          i = 1
          this%elem(i,j)%dxdeta = 0.5d0*(this%elem(i,j+1)%xc - this%elem(i,j-1)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(this%elem(i,j+1)%yc - this%elem(i,j-1)%yc)

          ! right side
          i = this%nelemi
          this%elem(i,j)%dxdeta = 0.5d0*(this%elem(i,j+1)%xc - this%elem(i,j-1)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(this%elem(i,j+1)%yc - this%elem(i,j-1)%yc)

        else if (j.eq.1) then

          ! southwest corner
          i = 1
          this%elem(i,j)%dxdeta = 0.5d0*(-3.0d0*this%elem(i,j)%xc + 4.0d0*this%elem(i,j+1)%xc - this%elem(i,j+2)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(-3.0d0*this%elem(i,j)%yc + 4.0d0*this%elem(i,j+1)%yc - this%elem(i,j+2)%yc)

          ! southeast corner
          i = this%nelemi
          this%elem(i,j)%dxdeta = 0.5d0*(-3.0d0*this%elem(i,j)%xc + 4.0d0*this%elem(i,j+1)%xc - this%elem(i,j+2)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(-3.0d0*this%elem(i,j)%yc + 4.0d0*this%elem(i,j+1)%yc - this%elem(i,j+2)%yc)

        else

          ! northwest corner
          i = 1
          this%elem(i,j)%dxdeta = 0.5d0*(3.0d0*this%elem(i,j)%xc - 4.0d0*this%elem(i,j-1)%xc + this%elem(i,j-2)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(3.0d0*this%elem(i,j)%yc - 4.0d0*this%elem(i,j-1)%yc + this%elem(i,j-2)%yc)

          ! northeast corner
          i = this%nelemi
          this%elem(i,j)%dxdeta = 0.5d0*(3.0d0*this%elem(i,j)%xc - 4.0d0*this%elem(i,j-1)%xc + this%elem(i,j-2)%xc)
          this%elem(i,j)%dydeta = 0.5d0*(3.0d0*this%elem(i,j)%yc - 4.0d0*this%elem(i,j-1)%yc + this%elem(i,j-2)%yc)

        end if


      end do

      ! Computing metrics on the horizontal boundaries
      ! PARALLELIZE
      do i=2,this%nelemi-1

        ! top
        j = this%nelemj
        this%elem(i,j)%dxdxi  = 0.5d0*(this%elem(i+1,j)%xc - this%elem(i-1,j)%xc)
        this%elem(i,j)%dydxi  = 0.5d0*(this%elem(i+1,j)%yc - this%elem(i-1,j)%yc)
        this%elem(i,j)%dxdeta = 0.5d0*(3.0d0*this%elem(i,j)%xc - 4.0d0*this%elem(i,j-1)%xc + this%elem(i,j-2)%xc)
        this%elem(i,j)%dydeta = 0.5d0*(3.0d0*this%elem(i,j)%yc - 4.0d0*this%elem(i,j-1)%yc + this%elem(i,j-2)%yc)

        ! bottom
        j = 1
        this%elem(i,j)%dxdxi  = 0.5d0*(this%elem(i+1,j)%xc - this%elem(i-1,j)%xc)
        this%elem(i,j)%dydxi  = 0.5d0*(this%elem(i+1,j)%yc - this%elem(i-1,j)%yc)
        this%elem(i,j)%dxdeta = 0.5d0*(-3.0d0*this%elem(i,j)%xc + 4.0d0*this%elem(i,j+1)%xc - this%elem(i,j+2)%xc)
        this%elem(i,j)%dydeta = 0.5d0*(-3.0d0*this%elem(i,j)%yc + 4.0d0*this%elem(i,j+1)%yc - this%elem(i,j+2)%yc)

      end do

      ! Computing Jacobian determinant for all interior elements
      do j=1,this%nelemj
        do i=1,this%nelemi
          this%elem(i,j)%detJ = this%elem(i,j)%dxdxi*this%elem(i,j)%dydeta - this%elem(i,j)%dxdeta*this%elem(i,j)%dydxi
        end do
      end do

      ! Finding lengths of all the edges
      ! Horizontal edges
      do j=1,this%jmax
        do i=1,this%nelemi

          ! Getting node locations
          x1 = this%x(i,j)
          y1 = this%y(i,j)
          x2 = this%x(i+1,j)
          y2 = this%y(i+1,j)

          ! Computing distance between the two nodes
          this%edges_h(i,j)%length = sqrt((x2-x1)**2 + (y2-y1)**2)

          ! Computing the midpoint of the edge
          this%edges_h(i,j)%xm = (x1+x2)/2.0d0
          this%edges_h(i,j)%ym = (y1+y2)/2.0d0

        end do
      end do

      ! Vertical edges
      do j=1,this%nelemj
        do i=1,this%imax

          ! Getting node locations
          x1 = this%x(i,j)
          y1 = this%y(i,j)
          x2 = this%x(i,j+1)
          y2 = this%y(i,j+1)

          ! Computing distance between the two nodes
          this%edges_v(i,j)%length = sqrt((x2-x1)**2 + (y2-y1)**2)

          ! Computing the midpoint of the edge
          this%edges_v(i,j)%xm = (x1+x2)/2.0d0
          this%edges_v(i,j)%ym = (y1+y2)/2.0d0

        end do
      end do

      print *, "Done preprocessing mesh."

    end subroutine preprocess

end module mesh_class
