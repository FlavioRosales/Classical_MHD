module hydro_diagnostics_lib
  use hydro_base_lib
  implicit none

  type, extends(hydro_base) :: hydro_diagnostics
  contains

    procedure :: rhomax
    procedure :: mass, mass_center
    procedure :: px_lineal, py_lineal, pz_lineal
    procedure :: Lx_angular, Ly_angular, Lz_angular
    procedure :: K_energy, W_energy, U_energy

  end type hydro_diagnostics

contains

  !=============================================================================!
  real(kind=8) function rhomax(this)
    implicit none
    class(hydro_diagnostics), intent(in out) :: this

    call MPI_ALLREDUCE(maxval(this%rho), rhomax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

  end function rhomax
  !=============================================================================!
  real(kind=8) function mass(this, fourier_space)
    use integrals
    use fft_lib
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    logical, intent(in) :: fourier_space
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: f

    if(fourier_space) then
      f = abs(fftshift(dcmplx(sqrt(this%rho))) * domain%dx * domain%dy * domain%dz)**2
      mass = trapezium(f, domain%dkx, domain%dky, domain%dkz)
    else
      mass = trapezium(this%rho, domain%dx, domain%dy, domain%dz)
    end if

  end function mass
  !============================================================================!
  function px_lineal(this) result(p)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    integer :: N
    real(kind=8), dimension(2) :: p
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: px_integrand

    N = domain%Nz/2

    px_integrand = dreal(transpose_y_to_x(dcmplx(this%rho * this%vx)))
    p(1) = trapezium(px_integrand(:,:,1:N), domain%dx, domain%dy, domain%dz)
    p(2) = trapezium(px_integrand(:,:,N+1:2*N), domain%dx, domain%dy, domain%dz)

  end function px_lineal
  !============================================================================!
  function py_lineal(this) result(p)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    integer :: N
    real(kind=8), dimension(2) :: p
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: py_integrand

    N = domain%Nz/2

    py_integrand = dreal(transpose_z_to_y(dcmplx(this%rho * this%vy)))
    p(1) = trapezium(py_integrand(:,:,1:N), domain%dx, domain%dy, domain%dz)
    p(2) = trapezium(py_integrand(:,:,N+1:2*N), domain%dx, domain%dy, domain%dz)

  end function py_lineal
  !============================================================================!
  function pz_lineal(this) result(p)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    integer :: N
    real(kind=8), dimension(2) :: p
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: pz_integrand

    N = domain%Nz/2

    pz_integrand = this%rho * this%vx
    p(1) = trapezium(pz_integrand(:,:,1:N), domain%dx, domain%dz, domain%dy)
    p(2) = trapezium(pz_integrand(:,:,N+1:2*N), domain%dx, domain%dz, domain%dy)

  end function pz_lineal
  !============================================================================!
  real(kind=8) function Lx_angular(this)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this

    Lx_angular = trapezium(this%rho * (domain%y * this%vz - domain%z * this%vy), &
    domain%dx, domain%dy, domain%dz)

  end function Lx_angular
  !============================================================================!
  real(kind=8) function Ly_angular(this)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this

    Ly_angular = trapezium(this%rho * (domain%z * this%vx - domain%x * this%vz), &
    domain%dx, domain%dy, domain%dz)

  end function Ly_angular
  !============================================================================!
  real(kind=8) function Lz_angular(this)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this

    Lz_angular = trapezium(this%rho * (domain%x * this%vy - domain%y * this%vx), &
    domain%dx, domain%dy, domain%dz)

  end function Lz_angular
  !============================================================================!
  real(kind=8) function K_energy(this, fourier_space)
    use integrals
    use fft_lib, only: fftshift
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    logical, intent(in) :: fourier_space

    if(fourier_space) then
      K_energy = trapezium(abs(fftshift(dcmplx(dsqrt(0.50d0 * this%rho * &
      (this%vx**2 + this%vy**2 + this%vz**2)) &
      * domain%dx * domain%dy * domain%dz)))**2, domain%dkx, domain%dky, domain%dkz)
    else
      K_energy = trapezium(0.50d0 * this%rho * (this%vx**2 + this%vy**2 + this%vz**2), domain%dx, domain%dy, domain%dz)
    end if

  end function K_energy
  !============================================================================!
  real(kind=8) function W_energy(this, fourier_space)
    use integrals
    use fft_lib, only: fftshift
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    logical, intent(in) :: fourier_space

    if(fourier_space) then
      W_energy = trapezium(- 0.50d0 * abs(fftshift(dcmplx(dsqrt(abs(this%V * this%rho)))) &
      * domain%dx * domain%dy * domain%dz)**2, domain%dkx, domain%dky, domain%dkz)
    else
      W_energy = trapezium(+ 0.50d0 * this%V * this%rho, domain%dx, domain%dy, domain%dz)
    end if

  end function W_energy
  !============================================================================!
  real(kind=8) function U_energy(this, fourier_space)
    use integrals
    use fft_lib, only: fftshift
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    logical, intent(in) :: fourier_space

    if(fourier_space) then
      U_energy = trapezium(abs(fftshift(dcmplx(dsqrt(this%rho *this%e) &
      * domain%dx * domain%dy * domain%dz)))**2, domain%dkx, domain%dky, domain%dkz)
    else
      U_energy = trapezium(this%rho * this%e, domain%dx, domain%dy, domain%dz)
    end if

  end function U_energy
  !============================================================================!
  function mass_center(this) result(xcm)
    use integrals
    implicit none
    class(hydro_diagnostics), intent(in out) :: this
    real(kind=8), dimension(3) :: xcm

    xcm(1) = trapezium(domain%x * this%rho, domain%dx, domain%dy, domain%dz)
    xcm(2) = trapezium(domain%y * this%rho, domain%dx, domain%dy, domain%dz)
    xcm(3) = trapezium(domain%z * this%rho, domain%dx, domain%dy, domain%dz)

  end function mass_center
  !============================================================================!

end module hydro_diagnostics_lib
