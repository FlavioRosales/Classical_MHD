module mesh_lib
  use mpi_lib
  implicit none

  type, public :: mesh

    integer :: Nx, Ny, Nz, Nt
    real(kind=8) :: xmin, xmax, dx, dkx
    real(kind=8) :: ymin, ymax, dy, dky
    real(kind=8) :: zmin, zmax, dz, dkz
    real(kind=8) :: tmin, tmax, dt, CFL, t, t_p
    real(kind=8), allocatable, dimension(:,:,:) :: x, y, z, r
    real(kind=8), allocatable, dimension(:,:,:) :: kx, ky, kz, k_square
    real(kind=8), allocatable, dimension(:,:,:) :: px, py, pz, p_square
    logical :: save_x, save_y, save_z, save_xy, save_xz, save_yz
    integer :: i_save, j_save, k_save

  contains

    procedure :: save_indexes
    procedure, private :: memory
    procedure, private :: create

  end type mesh

  interface mesh
    module procedure :: mesh_constructor
  end interface mesh

  type(mesh) :: domain

contains
  !============================================================================!
  function mesh_constructor(&
    xmin, xmax, Nx, ymin, ymax, Ny, zmin, zmax, Nz, tmin, tmax, CFL) result(this)
    implicit none
    !
    ! input
    !
    integer, intent(in) :: Nx, Ny, Nz
    real(kind=8), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax, tmin, tmax, CFL
    !
    ! output
    !
    type(mesh) :: this

    if(mod(Nx,row_nproc).ne.0) stop "The number of processors must divide the number of cells."
    if(mod(Ny,col_nproc).ne.0) stop "The number of processors must divide the number of cells."

    this%Nx = Nx/row_nproc; this%Ny = Ny/col_nproc; this%Nz = Nz

    this%xmin = xmin; this%xmax = xmax; this%dx = (xmax - xmin) / dble(Nx-1)
    this%ymin = ymin; this%ymax = ymax; this%dy = (ymax - ymin) / dble(Ny-1)
    this%zmin = zmin; this%zmax = zmax; this%dz = (zmax - zmin) / dble(Nz-1)

    this%tmin = tmin; this%tmax = tmax; this%CFL = CFL; this%t = tmin
    this%dt = CFL * min(this%dx, this%dy, this%dz)
    this%Nt = int((this%tmax - this%tmin) / this%dt)
    if(mod(this%Nt,2).ne.0) this%Nt = this%Nt + 1

    call this%memory
    call this%create

  end function mesh_constructor
  !============================================================================!
  subroutine save_indexes(this, x_save, y_save, z_save)
    implicit none
    !
    ! input
    !
    class(mesh), intent(in out) :: this
    real(kind=8), intent(in) :: x_save, y_save, z_save

    this%save_x = .false.; this%save_y = .false.; this%save_z = .false.

    if(this%x(1,1,1) <= x_save - this%dx .and. x_save - this%dx <= this%x(this%Nx,1,1)) this%save_x = .true.
    if(this%y(1,1,1) <= y_save - this%dy .and. y_save - this%dy <= this%y(1,this%Ny,1)) this%save_y = .true.
    if(this%z(1,1,1) <= z_save - this%dz .and. z_save - this%dz <= this%z(1,1,this%Nz)) this%save_z = .true.

    if(this%save_x) this%i_save = int((x_save - this%x(1,1,1))/this%dx) + 1
    if(this%save_y) this%j_save = int((y_save - this%y(1,1,1))/this%dy) + 1
    if(this%save_z) this%k_save = int((z_save - this%z(1,1,1))/this%dz) + 1

    this%save_xy = this%save_x .and. this%save_y
    this%save_xz = this%save_x .and. this%save_z
    this%save_yz = this%save_y .and. this%save_z

  end subroutine save_indexes
  !============================================================================!
  subroutine memory(this)
    implicit none
    !
    ! input
    !
    class(mesh), intent(in out) :: this

    allocate(this%x(this%Nx,this%Ny,this%Nz))
    allocate(this%y(this%Nx,this%Ny,this%Nz))
    allocate(this%z(this%Nx,this%Ny,this%Nz))
    allocate(this%r(this%Nx,this%Ny,this%Nz))
    allocate(this%kx(this%Nx,this%Ny,this%Nz))
    allocate(this%ky(this%Nx,this%Ny,this%Nz))
    allocate(this%kz(this%Nx,this%Ny,this%Nz))
    allocate(this%k_square(this%Nx,this%Ny,this%Nz))
    allocate(this%px(this%Nx,this%Ny,this%Nz))
    allocate(this%py(this%Nx,this%Ny,this%Nz))
    allocate(this%pz(this%Nx,this%Ny,this%Nz))
    allocate(this%p_square(this%Nx,this%Ny,this%Nz))

  end subroutine memory
  !============================================================================!
  subroutine create(this)
    !
    ! input
    !
    class(mesh), intent(in out) :: this
    !
    ! internal
    !
    integer :: i, j, k
    real(kind=8) :: xmin, ymin, zmin
    real(kind=8) :: dkx, dky, dkz
    real(kind=8), parameter :: pi_2 = 2.0d0 * acos(-1.0d0)

    xmin = this%xmin + dble(coords(1) * (this%Nx)) * this%dx
    ymin = this%ymin + dble(coords(2) * (this%Ny)) * this%dy
    zmin = this%zmin

    do i=1, this%Nx
      do j=1, this%Ny
        do k=1, this%Nz

          this%x(i,j,k) = xmin + dble(i-1) * this%dx
          this%y(i,j,k) = ymin + dble(j-1) * this%dy
          this%z(i,j,k) = zmin + dble(k-1) * this%dz

        end do
      end do
    end do

    this%r = sqrt(this%x**2 + this%y**2 + this%z**2)

    dkx = pi_2 / (dble(this%Nz) * this%dx)
    dky = pi_2 / (dble(this%Nz) * this%dy)
    dkz = pi_2 / (dble(this%Nz) * this%dz)

    this%dkx = 1.0d0 / (dble(this%Nz) * this%dx)
    this%dky = 1.0d0 / (dble(this%Nz) * this%dy)
    this%dkz = 1.0d0 / (dble(this%Nz) * this%dz)

    do i=1, this%Nx
      do j=1, this%Ny
        do k=1, this%Nz

          this%px(i,j,k) = dble(k-this%Nz/2-1) * dkx
          this%py(i,j,k) = dble(k-this%Nz/2-1) * dky
          this%pz(i,j,k) = dble(k-this%Nz/2-1) * dkz

        end do
      end do
    end do

    this%py = dreal(transpose_y_to_x(dcmplx(this%py)))
    this%pz = dreal(transpose_z_to_y(dcmplx(this%pz)))

    this%p_square = this%px**2 + this%py**2 + this%pz**2

    do i=1, this%Nx
      do j=1, this%Ny
        do k=1, this%Nz

          this%kx(i,j,k) = dble(k-this%Nz/2-1) * dkx
          this%ky(i,j,k) = dble(k-this%Nz/2-1) * dky
          this%kz(i,j,k) = dble(k-this%Nz/2-1) * dkz

        end do
      end do
    end do

    this%k_square = this%kx
    this%kx(:,:,1:this%Nz/2) = this%k_square(:,:,this%Nz/2+1:this%Nz)
    this%kx(:,:,this%Nz/2+1:this%Nz) = this%k_square(:,:,1:this%Nz/2)

    this%k_square = this%ky
    this%ky(:,:,1:this%Nz/2) = this%k_square(:,:,this%Nz/2+1:this%Nz)
    this%ky(:,:,this%Nz/2+1:this%Nz) = this%k_square(:,:,1:this%Nz/2)

    this%k_square = this%kz
    this%kz(:,:,1:this%Nz/2) = this%k_square(:,:,this%Nz/2+1:this%Nz)
    this%kz(:,:,this%Nz/2+1:this%Nz) = this%k_square(:,:,1:this%Nz/2)

    this%ky = dreal(transpose_y_to_x(dcmplx(this%ky)))
    this%kz = dreal(transpose_z_to_y(dcmplx(this%kz)))

    this%k_square = this%kx**2 + this%ky**2 + this%kz**2

  end subroutine create
  !============================================================================!
end module mesh_lib
