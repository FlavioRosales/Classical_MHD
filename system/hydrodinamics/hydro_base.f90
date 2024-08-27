module hydro_base_lib
  use mesh_lib
  implicit none

  real(kind=8), parameter :: zero = 0.0d0
  real(kind=8), parameter ::  one = 1.0d0
  real(kind=8), parameter :: half = 0.50d0
  real(kind=8), parameter ::   pi = dacos(-1.0d0)
  complex(kind=8), parameter, private :: i_unreal = (0.0d0, 1.0d0)

  real(kind=8), private :: rhoL, vxL, vyL, vzL, eL, pL, BxL, ByL, BzL, BL, vL, BvL 
  real(kind=8), private :: CaxL, CayL, CazL 
  real(kind=8), private :: CfxL, CfyL, CfzL 
  real(kind=8), private :: CsxL, CsyL, CszL  
  real(kind=8), private :: rhoR, vxR, vyR, vzR, eR, pR, BxR, ByR, BzR, BR, vR, BvR
  real(kind=8), private :: CaxR, CayR, CazR 
  real(kind=8), private :: CfxR, CfyR, CfzR 
  real(kind=8), private :: CsxR, CsyR, CszR 

  type, public :: hydro_base

    real(kind=8) :: gamma, floor, K_poly

    real(kind=8), allocatable, dimension(:,:,:) :: rho, vx, vy, vz, e, p, V
    real(kind=8), allocatable, dimension(:,:,:,:) :: u, u_p, Fx, Fy, Fz, S

    character(len=:), allocatable :: EoS, flux_formula, boundary

  contains

    procedure :: hydro_memory
    procedure :: write_check_point
    procedure :: read_check_point
    procedure :: primitive_variables
    procedure :: conservative_variables
    procedure :: numerical_flux
    procedure :: numerical_source
    procedure :: numerical_flux_derivatives
    procedure :: reconstructor, MINMOD_periodic, MINMOD_isolated
    procedure :: HLLE_x, HLLE_y, HLLE_z
    procedure :: MARQUINA_x, MARQUINA_y, MARQUINA_z
    procedure :: rhs_hydro
    procedure :: evolve
    procedure :: boundary_isolated

  end type hydro_base

contains
  !============================================================================!
  subroutine hydro_memory(this)
    implicit none
    class(hydro_base), intent(in out) :: this

    allocate(this% rho(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  vx(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  vy(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  vz(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  Bx(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  By(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  Bz(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   e(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   p(domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   V(domain%Nx, domain%Ny, domain%Nz))

    allocate(this%    u(8, domain%Nx, domain%Ny, domain%Nz))
    allocate(this%  u_p(8, domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   Fx(8, domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   Fy(8, domain%Nx, domain%Ny, domain%Nz))
    allocate(this%   Fz(8, domain%Nx, domain%Ny, domain%Nz))
    allocate(this%    S(8, domain%Nx, domain%Ny, domain%Nz))

  end subroutine hydro_memory
  !============================================================================!
  subroutine write_check_point(this, name, time)
    use strings_lib
    implicit none
    class(hydro_base), intent(in out) :: this
    character(len=*), intent(in) :: name
    real(kind=8), intent(in) :: time

    call system('mkdir -p check_point')

    open(123, file = 'check_point/' // name // char_integer(rank+1), action = 'write')
      write(123,*) time
      write(123,*) this%rho, this%vx, this%vy, this%vz, this%e
    close(123)

  end subroutine write_check_point
  !============================================================================!
  subroutine read_check_point(this, name, time)
    use strings_lib
    implicit none
    class(hydro_base), intent(in out) :: this
    character(len=*), intent(in) :: name
    real(kind=8), intent(in out) :: time


    open(123, file = 'check_point/'// name // char_integer(rank+1), action = 'read')
      read(123,*) time
      read(123,*) this%rho, this%vx, this%vy, this%vz, this%e
    close(123)

  end subroutine read_check_point
  !============================================================================!
  subroutine evolve(this)
    implicit none
    class(hydro_base), intent(in out) :: this

    integer :: rk

    this%u_p = this%u

    do rk=1, 3

      if(rk.eq.1) then
        this%u= this%u_p + this%rhs_hydro() * domain%dt
      else if(rk.eq.2) then
        this%u= 0.750d0 * this%u_p &
        + 0.250d0 * (this%u+ this%rhs_hydro() * domain%dt)
      else
        this%u= this%u_p/3.0d0  + 2.0d0 * (this%u + this%rhs_hydro() * domain%dt) / 3.0d0
      end if

    end do

  end subroutine evolve
  !============================================================================!
  subroutine primitive_variables(this)
    implicit none
    class(hydro_base), intent(in out) :: this

    this%u(1,:,:,:) = max(this%u(1,:,:,:), this%floor)

    this%rho = this%u(1,:,:,:)
    this% vx = this%u(2,:,:,:) / this%u(1,:,:,:)
    this% vy = this%u(3,:,:,:) / this%u(1,:,:,:)
    this% vz = this%u(4,:,:,:) / this%u(1,:,:,:)
    this% Bx = this%u(6,:,:,:) 
    this% By = this%u(7,:,:,:) 
    this% Bz = this%u(8,:,:,:) 
    this%  e = this%u(5,:,:,:) / this%u(1,:,:,:) - half * (this%vx**2 + this%vy**2 + this%vz**2) &
               - half*(this%Bx**2 + this%By**2 + this%Bz**2)
    if(this%EoS.eq.'ideal gas') then
      this%p = (this%gamma - one) * this%e * this%rho
    else if(this%EoS.eq.'poly') then
      this%p = this%K_poly * this%rho**this%gamma
      this%e = this%p / this%rho / (this%gamma - one)
    else if(this%EoS.eq.'dust') then
      this%p = zero
      this%e = zero
    end if

  end subroutine primitive_variables
  !============================================================================!
  subroutine conservative_variables(this)
    implicit none 
    class(hydro_base), intent(in out) :: this 

    this%u(1,:,:,:) = this%rho
    this%u(2,:,:,:) = this%rho * this%vx
    this%u(3,:,:,:) = this%rho * this%vy
    this%u(4,:,:,:) = this%rho * this%vz
    this%u(5,:,:,:) = this%e*this%rho + half*this%rho*(this%vx**2 + this%vy**2 + this%vz**2) &
                      + half*(this%Bx**2 + this%By**2 + this%Bz**2)
    this%u(6,:,:,:) = this%Bx
    this%u(7,:,:,:) = this%By
    this%u(8,:,:,:) = this%Bz

  end subroutine conservative_variables
  !============================================================================!
  function rhs_hydro(this) result(rhs)
    use mpi_lib, only: transpose_y_to_x_real, transpose_z_to_y_real
    implicit none
    class(hydro_base), intent(in out) :: this
    real(kind=8), dimension(8, domain%Nx, domain%Ny, domain%Nz) :: rhs

    call this%primitive_variables
    call this%numerical_source
    call this%numerical_flux
    call this%numerical_flux_derivatives

    rhs = this%S

  end function rhs_hydro
  !============================================================================!
  subroutine numerical_flux(this)
    implicit none
    class(hydro_base), intent(in out) :: this

    if(this%flux_formula.eq.'HLLE') then
      call this%HLLE_x
      call this%HLLE_y
      call this%HLLE_z
    else if(this%flux_formula.eq.'MARQUINA') then
      call this%MARQUINA_x
      call this%MARQUINA_y
      call this%MARQUINA_z
    else
      stop "The flux formula can only be 'HLLE' or 'MARQUINA'."
    end if

  end subroutine numerical_flux
  !============================================================================!
  subroutine numerical_source(this)
    use finite_differences, only: &
    dx_isolated => first_derivative_x_2_real, &
    dy_isolated => first_derivative_y_2_real, &
    dz_isolated => first_derivative_z_2_real
    use periodic_finite_differences, only: &
    dx_periodic => first_derivative_x_2_real, &
    dy_periodic => first_derivative_y_2_real, &
    dz_periodic => first_derivative_z_2_real
    implicit none
    class(hydro_base), intent(in out) :: this
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: dxV, dyV, dzV

    if(this%boundary.eq.'isolated') then
      dxV = dx_isolated(this%V, domain%dx)
      dyV = dy_isolated(this%V, domain%dy)
      dzV = dz_isolated(this%V, domain%dz)
    else if(this%boundary.eq.'periodic') then
      dxV = dx_periodic(this%V, domain%dx)
      dyV = dy_periodic(this%V, domain%dy)
      dzV = dz_periodic(this%V, domain%dz)
    end if

    this%S(1,:,:,:) = zero
    this%S(2,:,:,:) = - this%rho * dxV
    this%S(3,:,:,:) = - this%rho * dyV
    this%S(4,:,:,:) = - this%rho * dzV
    this%S(5,:,:,:) = - this%rho * (this%vx * dxV + this%vy * dyV + this%vz * dzV)

  end subroutine numerical_source
  !============================================================================!
  subroutine numerical_flux_derivatives(this)
    implicit none
    class(hydro_base), intent(in out) :: this
    real(kind=8) :: &
    dFx(5, domain%Nx, domain%Ny, domain%Nz), &
    dFy(5, domain%Nx, domain%Ny, domain%Nz), &
    dFz(5, domain%Nx, domain%Ny, domain%Nz)

    integer :: k

    dFx(:,:,:,1) = (this%Fx(:,:,:,1) - this%Fx(:,:,:,domain%Nz)) / domain%dx
    do k=2, domain%Nz
      dFx(:,:,:,k) = (this%Fx(:,:,:,k) - this%Fx(:,:,:,k-1)) / domain%dx
    end do
    do k=1, 5
      dFx(k,:,:,:) = transpose_y_to_x(dFx(k,:,:,:))
    end do
    !this%S = this%S - dF

    dFy(:,:,:,1) = (this%Fy(:,:,:,1) - this%Fy(:,:,:,domain%Nz)) / domain%dy
    do k=2, domain%Nz
      dFy(:,:,:,k) = (this%Fy(:,:,:,k) - this%Fy(:,:,:,k-1)) / domain%dy
    end do
    do k=1, 5
      dFy(k,:,:,:) = transpose_z_to_y(dFy(k,:,:,:))
    end do
    !this%S = this%S - dF

    dFz(:,:,:,1) = (this%Fz(:,:,:,1) - this%Fz(:,:,:,domain%Nz)) / domain%dz
    do k=2, domain%Nz
      dFz(:,:,:,k) = (this%Fz(:,:,:,k) - this%Fz(:,:,:,k-1)) / domain%dz
    end do
    this%S = this%S - dFx - dFy - dFz

  end subroutine numerical_flux_derivatives
  !============================================================================!
  subroutine HLLE_x(this)
    use mpi_lib, only: transpose_y_to_x
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k
    real(kind=8), dimension(8) :: lambdaL, lambdaR, FL, FR, uL, uR

    this%rho = transpose_y_to_x(this%rho)
    this% vx = transpose_y_to_x(this% vx)
    this% vy = transpose_y_to_x(this% vy)
    this% vz = transpose_y_to_x(this% vz)
    this% Bx = transpose_y_to_x(this% Bx)
    this% By = transpose_y_to_x(this% By)
    this% Bz = transpose_y_to_x(this% Bz)
    this%  e = transpose_y_to_x(this%  e)

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          lambdaL = [vxL - CfxL, &
          vxL - CaxL, &
          vxL - CsxL, & 
          vxL, & 
          vxL, & 
          vxL + CsxL, & 
          vxL + CaxL, &
          vxL + CfxL]

          lambdaR = [vxR - CfxR, &
          vxR - CaxR, &
          vxR - CsxR, & 
          vxR, & 
          vxR, & 
          vxR + CsxR, & 
          vxR + CaxR, &
          vxR + CfxR]

          uL = [rhoL, &
          rhoL*vxL, &
          rhoL*vyL, &
          rhoL*vzL, &
          eL*rhoL + half*rhoL*vL + half*BL, &
          BxL, & 
          ByL, &
          BzL]

          uR = [rhoR, &
          rhoR*vxR, &
          rhoR*vyR, &
          rhoR*vzR, &
          eR*rhoR + half*rhoR*vR + half*BR, &
          BxR, & 
          ByR, &
          BzR]


          SL = min(zero, minval(lambdaL), minval(lambdaR))
          SR = max(zero, maxval(lambdaL), maxval(lambdaR))

          FL = [uL(2), &
          uL(2) * vxL + pL + half*BL - BxL**2, &
          uL(2) * vyL - BxL*ByL, & 
          uL(2) * vzL - BxL*BzL, & 
          vxL * (uL(5) + pL + half*BL) - BvL*BxL, &
          zero, & 
          ByL*vxL - BxL*vyL, &
          BzL*vxL - BxL*vzL]

          FR = [uR(2), &
          uR(2) * vxR + pR + half*BR - BxR**2, &
          uR(2) * vyR - BxR*ByR, & 
          uR(2) * vzR - BxR*BzR, & 
          vxR * (uR(5) + pR + half*BR) - BvR*BxR, &
          zero, & 
          ByR*vxR - BxR*vyR, &
          BzR*vxR - BxR*vzR]

          if(SL.ne.SR) then
            this%Fx(:,i,j,k) = (SR * FL - SL * FR + SL * SR * (uR - uL)) / (SR - SL)
          else
            this%Fx(:,i,j,k) = half * (FL + FR - max(abs(SL), abs(SR)) * (uR - uL))
          end if

        end do
      end do
    end do

    this%rho = transpose_y_to_x(this%rho)
    this% vx = transpose_y_to_x(this% vx)
    this% vy = transpose_y_to_x(this% vy)
    this% vz = transpose_y_to_x(this% vz)
    this% Bx = transpose_y_to_x(this% Bx)
    this% By = transpose_y_to_x(this% By)
    this% Bz = transpose_y_to_x(this% Bz)
    this%  e = transpose_y_to_x(this%  e)

  end subroutine HLLE_x
  !============================================================================!
  subroutine HLLE_y(this)
    use mpi_lib, only: transpose_z_to_y
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k
    real(kind=8), dimension(8) :: lambdaL, lambdaR, FL, FR, uL, uR

    this%rho = transpose_z_to_y(this%rho)
    this% vx = transpose_z_to_y(this% vx)
    this% vy = transpose_z_to_y(this% vy)
    this% vz = transpose_z_to_y(this% vz)
    this% Bx = transpose_z_to_y(this% Bx)
    this% By = transpose_z_to_y(this% By)
    this% Bz = transpose_z_to_y(this% Bz)
    this%  e = transpose_z_to_y(this%  e)

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          lambdaL = [vyL - CfyL, &
          vyL - CayL, &
          vyL - CsyL, & 
          vyL, & 
          vyL, & 
          vyL + CsyL, & 
          vyL + CayL, &
          vyL + CfyL]

          lambdaR = [vyR - CfyR, &
          vyR - CayR, &
          vyR - CsyR, & 
          vyR, & 
          vyR, & 
          vyR + CsyR, & 
          vyR + CayR, &
          vyR + CfyR]

          uL = [rhoL, &
          rhoL*vxL, &
          rhoL*vyL, &
          rhoL*vzL, &
          eL*rhoL + half*rhoL*vL + half*BL, &
          BxL, & 
          ByL, &
          BzL]

          uR = [rhoR, &
          rhoR*vxR, &
          rhoR*vyR, &
          rhoR*vzR, &
          eR*rhoR + half*rhoR*vR + half*BR, &
          BxR, & 
          ByR, &
          BzR]


          SL = min(zero, minval(lambdaL), minval(lambdaR))
          SR = max(zero, maxval(lambdaL), maxval(lambdaR))

          FL = [uL(3), &
          uL(3) * vxL - ByL*BxL, &
          uL(3) * vyL + pL + half*BL - ByL**2, & 
          uL(3) * vzL - ByL*BzL, & 
          vyL * (uL(5) + pL + half*BL) - BvL*ByL, &
          BxL*vyL - ByL*vxL, &
          zero, & 
          BzL*vyL - ByL*vzL ]

          FR = [uR(3), &
          uR(3) * vxR - ByR*BxR, &
          uR(3) * vyR + pR + half*BR - ByR**2, & 
          uR(3) * vzR - ByR*BzR, & 
          vyR * (uR(5) + pR + half*BR) - BvR*ByR, &
          BxR*vyR - ByR*vxR, &
          zero, & 
          BzR*vyR - ByR*vzR ]

          if(SL.ne.SR) then
            this%Fy(:,i,j,k) = (SR * FL - SL * FR + SL * SR * (uR - uL)) / (SR - SL)
          else
            this%Fy(:,i,j,k) = half * (FL + FR - max(abs(SL), abs(SR)) * (uR - uL))
          end if

        end do
      end do
    end do

    this%rho = transpose_z_to_y(this%rho)
    this% vx = transpose_z_to_y(this% vx)
    this% vy = transpose_z_to_y(this% vy)
    this% vz = transpose_z_to_y(this% vz)
    this% Bx = transpose_z_to_y(this% Bx)
    this% By = transpose_z_to_y(this% By)
    this% Bz = transpose_z_to_y(this% Bz)
    this%  e = transpose_z_to_y(this%  e)

  end subroutine HLLE_y
  !============================================================================!
  subroutine HLLE_z(this)
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k
    real(kind=8), dimension(5) :: lambdaL, lambdaR, FL, FR

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          lambdaL = [vzL - CfzL, &
          vzL - CazL, &
          vzL - CszL, & 
          vzL, & 
          vzL, & 
          vzL + CszL, & 
          vzL + CazL, &
          vzL + CfzL]

          lambdaR = [vzR - CfzR, &
          vzR - CazR, &
          vzR - CszR, & 
          vzR, & 
          vzR, & 
          vzR + CszR, & 
          vzR + CazR, &
          vzR + CfzR]

          uL = [rhoL, &
          rhoL*vxL, &
          rhoL*vyL, &
          rhoL*vzL, &
          eL*rhoL + half*rhoL*vL + half*BL, &
          BxL, & 
          ByL, &
          BzL]

          uR = [rhoR, &
          rhoR*vxR, &
          rhoR*vyR, &
          rhoR*vzR, &
          eR*rhoR + half*rhoR*vR + half*BR, &
          BxR, & 
          ByR, &
          BzR]


          SL = min(zero, minval(lambdaL), minval(lambdaR))
          SR = max(zero, maxval(lambdaL), maxval(lambdaR))

          FL = [uL(4), &
          uL(4) * vxL - BzL*BxL, &
          uL(4) * vyL - BzL*ByL, & 
          uL(4) * vzL + pL + half*BL - BzL**2, & 
          vzL * (uL(5) + pL + half*BL) - BvL*BzL, &
          BxL*vzL - BzL*vxL, &
          ByL*vzL - BzL*vyL, &
          zero &
          ]

          FR = [uR(4), &
          uR(4) * vxR - BzR*BxR, &
          uR(4) * vyR - BzR*ByR, & 
          uR(4) * vzR + pR + half*BR - BzR**2, & 
          vzR * (uR(5) + pR + half*BR) - BvR*BzR, &
          BxR*vzR - BzR*vxR, &
          ByR*vzR - BzR*vyR, &
          zero &
          ]


          if(SL.ne.SR) then
            this%Fz(:,i,j,k) = (SR * FL - SL * FR + SL * SR * (uR - uL)) / (SR - SL)
          else
            this%Fz(:,i,j,k) = half * (FL + FR - max(abs(SL), abs(SR)) * (uR - uL))
          end if

        end do
      end do
    end do

  end subroutine HLLE_z
  !============================================================================!
  subroutine MARQUINA_x(this)
    use mpi_lib, only: transpose_y_to_x_real
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k, l
    real(kind=8), dimension(5) :: lambdaL, lambdaR
    real(kind=8), dimension(5,5) :: XL, XR
    real(kind=8), dimension(5) :: phiL, phiR, phi_minus, phi_plus, WL, WR

    this%rho = transpose_y_to_x_real(this%rho)
    this% vx = transpose_y_to_x_real(this% vx)
    this% vy = transpose_y_to_x_real(this% vy)
    this% vz = transpose_y_to_x_real(this% vz)
    this% Bx = transpose_y_to_x_real(this% Bx)
    this% By = transpose_y_to_x_real(this% By)
    this% Bz = transpose_y_to_x_real(this% Bz)
    this%  e = transpose_y_to_x_real(this%  e)

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          WL = [half, this%gamma - one, zero, zero, half] * rhoL / this%gamma
          WR = [half, this%gamma - one, zero, zero, half] * rhoR / this%gamma

          lambdaL = [vxL - csL, vxL, vxL, vxL, vxL + csL]
          lambdaR = [vxR - csR, vxR, vxR, vxR, vxR + csR]

          XL(1,:) = [one, vxL - csL, vyL, vzL, HL - vxL * csL]
          XR(1,:) = [one, vxR - csR, vyR, vzR, HR - vxR * csR]

          XL(2,:) = [one, vxL, vyL, vzL, vL]
          XR(2,:) = [one, vxR, vyR, vzR, vR]

          XL(3,:) = [zero, zero, one, zero, vyL]
          XR(3,:) = [zero, zero, one, zero, vyR]

          XL(4,:) = [zero, zero, zero, one, vzL]
          XR(4,:) = [zero, zero, zero, one, vzR]

          XL(5,:) = [one, vxL + csL, vyL, vzL, HL + vxL * csL]
          XR(5,:) = [one, vxR + csR, vyR, vzR, HR + vxR * csR]

          phiL = [half * (vxL - csL), vxL * (this%gamma - one), zero, zero, half * (vxL + csL) ] * rhoL / this%gamma
          phiR = [half * (vxR - csR), vxR * (this%gamma - one), zero, zero, half * (vxR + csR) ] * rhoR / this%gamma

          do l=1, 5
            if(lambdaL(l) * lambdaR(l) <= zero) then
              phi_minus(l) = half * (phiR(l) - max(abs(lambdaL(l)), abs(lambdaR(l))) * WR(l))
              phi_plus(l)  = half * (phiL(l) + max(abs(lambdaL(l)), abs(lambdaR(l))) * WL(l))
            else
              if(lambdaL(l) > zero) then
                phi_minus(l) = zero
                phi_plus(l)  = phiL(l)
              else
                phi_minus(l) = phiR(l)
                phi_plus(l)  = zero
              end if
            end if
          end do

          this%Fx(:,i,j,k) = &
          + phi_plus(1) * XL(1,:) + phi_minus(1) * XR(1,:) &
          + phi_plus(2) * XL(2,:) + phi_minus(2) * XR(2,:) &
          + phi_plus(3) * XL(3,:) + phi_minus(3) * XR(3,:) &
          + phi_plus(4) * XL(4,:) + phi_minus(4) * XR(4,:) &
          + phi_plus(5) * XL(5,:) + phi_minus(5) * XR(5,:)

        end do
      end do
    end do

    this%rho = transpose_y_to_x_real(this%rho)
    this% vx = transpose_y_to_x_real(this% vx)
    this% vy = transpose_y_to_x_real(this% vy)
    this% vz = transpose_y_to_x_real(this% vz)
    this% Bx = transpose_y_to_x_real(this% Bx)
    this% By = transpose_y_to_x_real(this% By)
    this% Bz = transpose_y_to_x_real(this% Bz)
    this%  e = transpose_y_to_x_real(this%  e)

  end subroutine MARQUINA_x
  !============================================================================!
  subroutine MARQUINA_y(this)
    use mpi_lib, only: transpose_z_to_y_real
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k, l
    real(kind=8), dimension(5) :: lambdaL, lambdaR
    real(kind=8), dimension(5,5) :: XL, XR
    real(kind=8), dimension(5) :: phiL, phiR, phi_minus, phi_plus, WL, WR

    this%rho = transpose_z_to_y_real(this%rho)
    this% vx = transpose_z_to_y_real(this% vx)
    this% vy = transpose_z_to_y_real(this% vy)
    this% vz = transpose_z_to_y_real(this% vz)
    this% Bx = transpose_z_to_y_real(this% Bx)
    this% By = transpose_z_to_y_real(this% By)
    this% Bz = transpose_z_to_y_real(this% Bz)
    this%  e = transpose_z_to_y_real(this%  e)

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          WL = [half, this%gamma - one, zero, zero, half] * rhoL / this%gamma
          WR = [half, this%gamma - one, zero, zero, half] * rhoR / this%gamma

          lambdaL = [vyL - csL, vyL, vyL, vyL, vyL + csL]
          lambdaR = [vyR - csR, vyR, vyR, vyR, vyR + csR]

          XL(1,:) = [one, vxL, vyL - csL, vzL, HL - vyL * csL]
          XR(1,:) = [one, vxR, vyR - csR, vzR, HR - vyR * csR]

          XL(2,:) = [one, vxL, vyL, vzL, vL]
          XR(2,:) = [one, vxR, vyR, vzR, vR]

          XL(3,:) = [zero, one, zero, zero, vxL]
          XR(3,:) = [zero, one, zero, zero, vxR]

          XL(4,:) = [zero, zero, zero, one, vzL]
          XR(4,:) = [zero, zero, zero, one, vzR]

          XL(5,:) = [one, vxL, vyL + csL, vzL, HL + vyL * csL]
          XR(5,:) = [one, vxR, vyR + csR, vzR, HR + vyR * csR]

          phiL = [half * (vyL - csL), vyL * (this%gamma - one), zero, zero, half * (vyL + csL) ] * rhoL / this%gamma
          phiR = [half * (vyR - csR), vyR * (this%gamma - one), zero, zero, half * (vyR + csR) ] * rhoR / this%gamma

          do l=1, 5
            if(lambdaL(l) * lambdaR(l) < zero) then
              phi_minus(l) = half * (phiR(l) - max(abs(lambdaL(l)), abs(lambdaR(l))) * WR(l))
              phi_plus(l)  = half * (phiL(l) + max(abs(lambdaL(l)), abs(lambdaR(l))) * WL(l))
            else
              if(lambdaL(l) > zero) then
                phi_minus(l) = zero
                phi_plus(l)  = phiL(l)
              else
                phi_minus(l) = phiR(l)
                phi_plus(l)  = zero
              end if
            end if
          end do

          this%Fy(:,i,j,k) = &
          + phi_plus(1) * XL(1,:) + phi_minus(1) * XR(1,:) &
          + phi_plus(2) * XL(2,:) + phi_minus(2) * XR(2,:) &
          + phi_plus(3) * XL(3,:) + phi_minus(3) * XR(3,:) &
          + phi_plus(4) * XL(4,:) + phi_minus(4) * XR(4,:) &
          + phi_plus(5) * XL(5,:) + phi_minus(5) * XR(5,:)

        end do
      end do
    end do

    this%rho = transpose_z_to_y_real(this%rho)
    this% vx = transpose_z_to_y_real(this% vx)
    this% vy = transpose_z_to_y_real(this% vy)
    this% vz = transpose_z_to_y_real(this% vz)
    this% Bx = transpose_z_to_y_real(this% Bx)
    this% By = transpose_z_to_y_real(this% By)
    this% Bz = transpose_z_to_y_real(this% Bz)
    this%  e = transpose_z_to_y_real(this%  e)

  end subroutine MARQUINA_y
  !============================================================================!
  subroutine MARQUINA_z(this)
    implicit none
    class(hydro_base), intent(in out) :: this
    integer :: i, j, k, l
    real(kind=8), dimension(5) :: lambdaL, lambdaR
    real(kind=8), dimension(5,5) :: XL, XR
    real(kind=8), dimension(5) :: phiL, phiR, phi_minus, phi_plus, WL, WR

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          call this%reconstructor(i,j,k)

          WL = [half, this%gamma - one, zero, zero, half] * rhoL / this%gamma
          WR = [half, this%gamma - one, zero, zero, half] * rhoR / this%gamma

          lambdaL = [vzL - csL, vzL, vzL, vzL, vzL + csL]
          lambdaR = [vzR - csR, vzR, vzR, vzR, vzR + csR]

          XL(1,:) = [one, vxL, vyL, vzL - csL, HL - vzL * csL]
          XR(1,:) = [one, vxR, vyR, vzR - csR, HR - vzR * csR]

          XL(2,:) = [one, vxL, vyL, vzL, vL]
          XR(2,:) = [one, vxR, vyR, vzR, vR]

          XL(3,:) = [zero, one, zero, zero, vxL]
          XR(3,:) = [zero, one, zero, zero, vxR]

          XL(4,:) = [zero, zero, one, zero, vyL]
          XR(4,:) = [zero, zero, one, zero, vyR]

          XL(5,:) = [one, vxL, vyL, vzL + csL, HL + vzL * csL]
          XR(5,:) = [one, vxR, vyR, vzR + csR, HR + vzR * csR]

          phiL = [half * (vzL - csL), vzL * (this%gamma - one), zero, zero, half * (vzL + csL) ] * rhoL / this%gamma
          phiR = [half * (vzR - csR), vzR * (this%gamma - one), zero, zero, half * (vzR + csR) ] * rhoR / this%gamma

          do l=1, 5
            if(lambdaL(l) * lambdaR(l) < zero) then
              phi_minus(l) = half * (phiR(l) - max(abs(lambdaL(l)), abs(lambdaR(l))) * WR(l))
              phi_plus(l)  = half * (phiL(l) + max(abs(lambdaL(l)), abs(lambdaR(l))) * WL(l))
            else
              if(lambdaL(l) > zero) then
                phi_minus(l) = zero
                phi_plus(l)  = phiL(l)
              else
                phi_minus(l) = phiR(l)
                phi_plus(l)  = zero
              end if
            end if
          end do

          this%Fz(:,i,j,k) = &
          + phi_plus(1) * XL(1,:) + phi_minus(1) * XR(1,:) &
          + phi_plus(2) * XL(2,:) + phi_minus(2) * XR(2,:) &
          + phi_plus(3) * XL(3,:) + phi_minus(3) * XR(3,:) &
          + phi_plus(4) * XL(4,:) + phi_minus(4) * XR(4,:) &
          + phi_plus(5) * XL(5,:) + phi_minus(5) * XR(5,:)

        end do
      end do
    end do

  end subroutine MARQUINA_z
  !============================================================================!
  subroutine reconstructor(this,i,j,k)
    implicit none
    class(hydro_base), intent(in out) :: this
    real(kind=8) :: EpsL, EpsR, HtL, HtR
    integer, intent(in) :: i, j, k

    if(this%boundary.eq.'periodic') then
      call this%MINMOD_periodic(i,j,k)
    else if(this%boundary.eq.'isolated') then
      call this%MINMOD_isolated(i,j,k)
    end if
    ! density
    rhoL = max(rhoL, this%floor)
    rhoR = max(rhoR, this%floor)
    ! pressure and internal energy
    if(this%EoS.eq.'ideal gas') then
      eL = max(eL, zero)
      eR = max(eR, zero)
      pL = (this%gamma - one) * rhoL * eL
      pR = (this%gamma - one) * rhoR * eR
    else if(this%EoS.eq.'poly') then
      pL = this%K_poly * rhoL ** this%gamma
      pR = this%K_poly * rhoR ** this%gamma
      eL = pL / rhoL / (this%gamma - one)
      eR = pR / rhoR / (this%gamma - one)
		else if(this%EoS.eq.'dust') then
			pL = zero
			pR = zero
      eL = zero
      eR = zero
    end if

    !B^2
    BL = BxL**2 + ByL**2 + BzL**2
    BR = BxR**2 + ByR**2 + BzR**2
    !V^2
    vL = vxL**2 + vyL**2 + vzL**2
    vR = vxR**2 + vyR**2 + vzR**2
    !B dot v
    BvL = BxL*vxL + ByL*vyL + BzL*vzL
    BvR = BxR*vxR + ByR*vyR + BzR*vzR
    !Alfven velocities
    CaxL = dsqrt((BxL)**2 / rhoL)
    CaxR = dsqrt((BxR)**2 / rhoR)
    CayL = dsqrt((ByL)**2 / rhoL)
    CayR = dsqrt((ByR)**2 / rhoR)
    CazL = dsqrt((BzL)**2 / rhoL)
    CazR = dsqrt((BzR)**2 / rhoR)
    ! Magnetosonic and acoustic velocities
    HtL = (this%gamma*pL + BL)/rhoL
    HtL = (this%gamma*pR + BR)/rhoR

    EpsL = HtL**2  &
           - 4.0d0*(this%gamma*pL*BxL**2)/rhoL**2
    EpsR = HtR**2  &
           - 4.0d0*(this%gamma*pR*BxR**2)/rhoR**2

    CfxL = dsqrt(half*(HtL + dsqrt(EpsL)))
    CSxL = dsqrt(half*(HtL - dsqrt(EpsL)))

    EpsL = HtL**2  &
    - 4.0d0*(this%gamma*pL*ByL**2)/rhoL**2
    EpsR = HtR**2  &
    - 4.0d0*(this%gamma*pR*ByR**2)/rhoR**2

    CfyL = dsqrt(half*(HtL + dsqrt(EpsL)))
    CSyL = dsqrt(half*(HtL - dsqrt(EpsL)))

    EpsL = HtL**2  &
    - 4.0d0*(this%gamma*pL*BzL**2)/rhoL**2
    EpsR = HtR**2  &
    - 4.0d0*(this%gamma*pR*BzR**2)/rhoR**2

    CfzL = dsqrt(half*(HtL + dsqrt(EpsL)))
    CSzL = dsqrt(half*(HtL - dsqrt(EpsL)))



  end subroutine reconstructor
  !============================================================================!
  subroutine MINMOD_periodic(this,i,j,k)
    implicit none
    class(hydro_base), intent(in out) :: this
    integer, intent(in) ::  i, j, k
    real(kind=8) :: delta_minus, delta_plus

    if(k>1 .and. k<domain%Nz-1) then

      ! density
      delta_minus = this%rho(i,j,k) - this%rho(i,j,k-1)
      delta_plus  = this%rho(i,j,k+1) - this%rho(i,j,k)
      rhoL = this%rho(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%rho(i,j,k+1) - this%rho(i,j,k)
      delta_plus  = this%rho(i,j,k+2) - this%rho(i,j,k+1)
      rhoR = this%rho(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! x-axis velocity
      delta_minus = this%vx(i,j,k) - this%vx(i,j,k-1)
      delta_plus  = this%vx(i,j,k+1) - this%vx(i,j,k)
      vxL = this%vx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vx(i,j,k+1) - this%vx(i,j,k)
      delta_plus  = this%vx(i,j,k+2) - this%vx(i,j,k+1)
      vxR = this%vx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis velocity
      delta_minus = this%vy(i,j,k) - this%vy(i,j,k-1)
      delta_plus  = this%vy(i,j,k+1) - this%vy(i,j,k)
      vyL = this%vy(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vy(i,j,k+1) - this%vy(i,j,k)
      delta_plus  = this%vy(i,j,k+2) - this%vy(i,j,k+1)
      vyR = this%vy(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis velocity
      delta_minus = this%vz(i,j,k) - this%vz(i,j,k-1)
      delta_plus  = this%vz(i,j,k+1) - this%vz(i,j,k)
      vzL = this%vz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vz(i,j,k+1) - this%vz(i,j,k)
      delta_plus  = this%vz(i,j,k+2) - this%vz(i,j,k+1)
      vzR = this%vz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! x-axis Magnetic field
      delta_minus = this%Bx(i,j,k) - this%Bx(i,j,k-1)
      delta_plus  = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      BxL = this%Bx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      delta_plus  = this%Bx(i,j,k+2) - this%Bx(i,j,k+1)
      BxR = this%Bx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis Magnetic field
      delta_minus = this%By(i,j,k) - this%By(i,j,k-1)
      delta_plus  = this%By(i,j,k+1) - this%By(i,j,k)
      ByL = this%By(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%By(i,j,k+1) - this%By(i,j,k)
      delta_plus  = this%By(i,j,k+2) - this%By(i,j,k+1)
      ByR = this%By(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis Magnetic field
      delta_minus = this%Bz(i,j,k) - this%Bz(i,j,k-1)
      delta_plus  = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      BzL = this%Bz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      delta_plus  = this%Bz(i,j,k+2) - this%Bz(i,j,k+1)
      BzR = this%Bz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! energy
      delta_minus = this%e(i,j,k) - this%e(i,j,k-1)
      delta_plus  = this%e(i,j,k+1) - this%e(i,j,k)
      eL = this%e(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%e(i,j,k+1) - this%e(i,j,k)
      delta_plus  = this%e(i,j,k+2) - this%e(i,j,k+1)
      eR = this%e(i,j,k+1) - half * delta(delta_minus, delta_plus)

    else if(k.eq.1) then

      ! density
      delta_minus = this%rho(i,j,k) - this%rho(i,j,domain%Nz)
      delta_plus  = this%rho(i,j,k+1) - this%rho(i,j,k)
      rhoL = max(this%rho(i,j,k) + half * delta(delta_minus, delta_plus), this%floor)

      delta_minus = this%rho(i,j,k+1) - this%rho(i,j,k)
      delta_plus  = this%rho(i,j,k+2) - this%rho(i,j,k+1)
      rhoR = max(this%rho(i,j,k+1) - half * delta(delta_minus, delta_plus), this%floor)
      ! x-axis velocity
      delta_minus = this%vx(i,j,k) - this%vx(i,j,domain%Nz)
      delta_plus  = this%vx(i,j,k+1) - this%vx(i,j,k)
      vxL = this%vx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vx(i,j,k+1) - this%vx(i,j,k)
      delta_plus  = this%vx(i,j,k+2) - this%vx(i,j,k+1)
      vxR = this%vx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis velocity
      delta_minus = this%vy(i,j,k) - this%vy(i,j,domain%Nz)
      delta_plus  = this%vy(i,j,k+1) - this%vy(i,j,k)
      vyL = this%vy(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vy(i,j,k+1) - this%vy(i,j,k)
      delta_plus  = this%vy(i,j,k+2) - this%vy(i,j,k+1)
      vyR = this%vy(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis velocity
      delta_minus = this%vz(i,j,k) - this%vz(i,j,domain%Nz)
      delta_plus  = this%vz(i,j,k+1) - this%vz(i,j,k)
      vzL = this%vz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vz(i,j,k+1) - this%vz(i,j,k)
      delta_plus  = this%vz(i,j,k+2) - this%vz(i,j,k+1)
      vzR = this%vz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! x-axis Magnetic field
      delta_minus = this%Bx(i,j,k) - this%Bx(i,j,domain%Nz)
      delta_plus  = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      BxL = this%Bx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      delta_plus  = this%Bx(i,j,k+2) - this%Bx(i,j,k+1)
      BxR = this%Bx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis Magnetic field
      delta_minus = this%By(i,j,k) - this%By(i,j,domain%Nz)
      delta_plus  = this%By(i,j,k+1) - this%By(i,j,k)
      ByL = this%By(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%By(i,j,k+1) - this%By(i,j,k)
      delta_plus  = this%By(i,j,k+2) - this%By(i,j,k+1)
      ByR = this%By(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis Magnetic field
      delta_minus = this%Bz(i,j,k) - this%Bz(i,j,domain%Nz)
      delta_plus  = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      BzL = this%Bz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      delta_plus  = this%Bz(i,j,k+2) - this%Bz(i,j,k+1)
      BzR = this%Bz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! energy
      delta_minus = this%e(i,j,k) - this%e(i,j,domain%Nz)
      delta_plus  = this%e(i,j,k+1) - this%e(i,j,k)
      eL = this%e(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%e(i,j,k+1) - this%e(i,j,k)
      delta_plus  = this%e(i,j,k+2) - this%e(i,j,k+1)
      eR = this%e(i,j,k+1) - half * delta(delta_minus, delta_plus)

    else if(k.eq.domain%Nz-1) then

      ! density
      delta_minus = this%rho(i,j,k) - this%rho(i,j,k-1)
      delta_plus  = this%rho(i,j,k+1) - this%rho(i,j,k)
      rhoL = max(this%rho(i,j,k) + half * delta(delta_minus, delta_plus), this%floor)

      delta_minus = this%rho(i,j,k+1) - this%rho(i,j,k)
      delta_plus  = this%rho(i,j,1) - this%rho(i,j,k+1)
      rhoR = max(this%rho(i,j,k+1) - half * delta(delta_minus, delta_plus), this%floor)
      ! x-axis velocity
      delta_minus = this%vx(i,j,k) - this%vx(i,j,k-1)
      delta_plus  = this%vx(i,j,k+1) - this%vx(i,j,k)
      vxL = this%vx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vx(i,j,k+1) - this%vx(i,j,k)
      delta_plus  = this%vx(i,j,1) - this%vx(i,j,k+1)
      vxR = this%vx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis velocity
      delta_minus = this%vy(i,j,k) - this%vy(i,j,k-1)
      delta_plus  = this%vy(i,j,k+1) - this%vy(i,j,k)
      vyL = this%vy(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vy(i,j,k+1) - this%vy(i,j,k)
      delta_plus  = this%vy(i,j,1) - this%vy(i,j,k+1)
      vyR = this%vy(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis velocity
      delta_minus = this%vz(i,j,k) - this%vz(i,j,k-1)
      delta_plus  = this%vz(i,j,k+1) - this%vz(i,j,k)
      vzL = this%vz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vz(i,j,k+1) - this%vz(i,j,k)
      delta_plus  = this%vz(i,j,1) - this%vz(i,j,k+1)
      vzR = this%vz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! x-axis Magnetic field
      delta_minus = this%Bx(i,j,k) - this%Bx(i,j,k-1)
      delta_plus  = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      BxL = this%Bx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bx(i,j,k+1) - this%Bx(i,j,k)
      delta_plus  = this%Bx(i,j,1) - this%Bx(i,j,k+1)
      BxR = this%Bx(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! y-axis Magnetic field
      delta_minus = this%By(i,j,k) - this%By(i,j,k-1)
      delta_plus  = this%By(i,j,k+1) - this%By(i,j,k)
      ByL = this%By(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%By(i,j,k+1) - this%By(i,j,k)
      delta_plus  = this%By(i,j,1) - this%By(i,j,k+1)
      ByR = this%By(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! z-axis Magnetic field
      delta_minus = this%Bz(i,j,k) - this%Bz(i,j,k-1)
      delta_plus  = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      BzL = this%Bz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bz(i,j,k+1) - this%Bz(i,j,k)
      delta_plus  = this%Bz(i,j,1) - this%Bz(i,j,k+1)
      BzR = this%Bz(i,j,k+1) - half * delta(delta_minus, delta_plus)
      ! energy
      delta_minus = this%e(i,j,k) - this%e(i,j,k-1)
      delta_plus  = this%e(i,j,k+1) - this%e(i,j,k)
      eL = this%e(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%e(i,j,k+1) - this%e(i,j,k)
      delta_plus  = this%e(i,j,1) - this%e(i,j,k+1)
      eR = this%e(i,j,k+1) - half * delta(delta_minus, delta_plus)

    else if(k.eq.domain%Nz) then

      ! density
      delta_minus = this%rho(i,j,k) - this%rho(i,j,k-1)
      delta_plus  = this%rho(i,j,1) - this%rho(i,j,k)
      rhoL = max(this%rho(i,j,k) + half * delta(delta_minus, delta_plus), this%floor)

      delta_minus = this%rho(i,j,1) - this%rho(i,j,k)
      delta_plus  = this%rho(i,j,2) - this%rho(i,j,1)
      rhoR = max(this%rho(i,j,1) - half * delta(delta_minus, delta_plus), this%floor)
      ! x-axis velocity
      delta_minus = this%vx(i,j,k) - this%vx(i,j,k-1)
      delta_plus  = this%vx(i,j,1) - this%vx(i,j,k)
      vxL = this%vx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vx(i,j,1) - this%vx(i,j,k)
      delta_plus  = this%vx(i,j,2) - this%vx(i,j,1)
      vxR = this%vx(i,j,2) - half * delta(delta_minus, delta_plus)
      ! y-axis velocity
      delta_minus = this%vy(i,j,k) - this%vy(i,j,k-1)
      delta_plus  = this%vy(i,j,1) - this%vy(i,j,k)
      vyL = this%vy(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vy(i,j,1) - this%vy(i,j,k)
      delta_plus  = this%vy(i,j,2) - this%vy(i,j,1)
      vyR = this%vy(i,j,1) - half * delta(delta_minus, delta_plus)
      ! z-axis velocity
      delta_minus = this%vz(i,j,k) - this%vz(i,j,k-1)
      delta_plus  = this%vz(i,j,1) - this%vz(i,j,k)
      vzL = this%vz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%vz(i,j,1) - this%vz(i,j,k)
      delta_plus  = this%vz(i,j,2) - this%vz(i,j,1)
      vzR = this%vz(i,j,2) - half * delta(delta_minus, delta_plus)
      ! x-axis Magnetic Field
      delta_minus = this%Bx(i,j,k) - this%Bx(i,j,k-1)
      delta_plus  = this%Bx(i,j,1) - this%Bx(i,j,k)
      BxL = this%Bx(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bx(i,j,1) - this%Bx(i,j,k)
      delta_plus  = this%Bx(i,j,2) - this%Bx(i,j,1)
      BxR = this%Bx(i,j,2) - half * delta(delta_minus, delta_plus)
      ! y-axis Magnetic Field
      delta_minus = this%By(i,j,k) - this%By(i,j,k-1)
      delta_plus  = this%By(i,j,1) - this%By(i,j,k)
      ByL = this%By(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%By(i,j,1) - this%By(i,j,k)
      delta_plus  = this%By(i,j,2) - this%By(i,j,1)
      ByR = this%By(i,j,1) - half * delta(delta_minus, delta_plus)
      ! z-axis Magnetic Field
      delta_minus = this%Bz(i,j,k) - this%Bz(i,j,k-1)
      delta_plus  = this%Bz(i,j,1) - this%Bz(i,j,k)
      BzL = this%Bz(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%Bz(i,j,1) - this%Bz(i,j,k)
      delta_plus  = this%Bz(i,j,2) - this%Bz(i,j,1)
      BzR = this%Bz(i,j,2) - half * delta(delta_minus, delta_plus)
      ! energy
      delta_minus = this%e(i,j,k) - this%e(i,j,k-1)
      delta_plus  = this%e(i,j,1) - this%e(i,j,k)
      eL = this%e(i,j,k) + half * delta(delta_minus, delta_plus)

      delta_minus = this%e(i,j,1) - this%e(i,j,k)
      delta_plus  = this%e(i,j,2) - this%e(i,j,1)
      eR = this%e(i,j,1) - half * delta(delta_minus, delta_plus)

    end if

  end subroutine MINMOD_periodic
  !============================================================================!
  subroutine MINMOD_isolated(this,i,j,k)
    implicit none
    class(hydro_base), intent(in out) :: this
    integer, intent(in) ::  i, j, k
    real(kind=8) :: delta_minus, delta_plus

    integer :: k0

    k0 = min(max(k, 2), domain%Nz-2)

    delta_minus = this%rho(i,j,k0) - this%rho(i,j,k0-1)
    delta_plus  = this%rho(i,j,k0+1) - this%rho(i,j,k0)
    rhoL = this%rho(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%rho(i,j,k0+1) - this%rho(i,j,k0)
    delta_plus  = this%rho(i,j,k0+2) - this%rho(i,j,k0+1)
    rhoR = this%rho(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! x-axis velocity
    delta_minus = this%vx(i,j,k0) - this%vx(i,j,k0-1)
    delta_plus  = this%vx(i,j,k0+1) - this%vx(i,j,k0)
    vxL = this%vx(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%vx(i,j,k0+1) - this%vx(i,j,k0)
    delta_plus  = this%vx(i,j,k0+2) - this%vx(i,j,k0+1)
    vxR = this%vx(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! y-axis velocity
    delta_minus = this%vy(i,j,k0) - this%vy(i,j,k0-1)
    delta_plus  = this%vy(i,j,k0+1) - this%vy(i,j,k0)
    vyL = this%vy(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%vy(i,j,k0+1) - this%vy(i,j,k0)
    delta_plus  = this%vy(i,j,k0+2) - this%vy(i,j,k0+1)
    vyR = this%vy(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! z-axis velocity
    delta_minus = this%vz(i,j,k0) - this%vz(i,j,k0-1)
    delta_plus  = this%vz(i,j,k0+1) - this%vz(i,j,k0)
    vzL = this%vz(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%vz(i,j,k0+1) - this%vz(i,j,k0)
    delta_plus  = this%vz(i,j,k0+2) - this%vz(i,j,k0+1)
    vzR = this%vz(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! x-axis Magnetic Field
    delta_minus = this%Bx(i,j,k0) - this%Bx(i,j,k0-1)
    delta_plus  = this%Bx(i,j,k0+1) - this%Bx(i,j,k0)
    BxL = this%Bx(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%Bx(i,j,k0+1) - this%Bx(i,j,k0)
    delta_plus  = this%Bx(i,j,k0+2) - this%Bx(i,j,k0+1)
    BxR = this%Bx(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! y-axis Magnetic Field
    delta_minus = this%By(i,j,k0) - this%By(i,j,k0-1)
    delta_plus  = this%By(i,j,k0+1) - this%By(i,j,k0)
    ByL = this%By(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%By(i,j,k0+1) - this%By(i,j,k0)
    delta_plus  = this%By(i,j,k0+2) - this%By(i,j,k0+1)
    ByR = this%By(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! z-axis Magnetic Field
    delta_minus = this%Bz(i,j,k0) - this%Bz(i,j,k0-1)
    delta_plus  = this%Bz(i,j,k0+1) - this%Bz(i,j,k0)
    BzL = this%Bz(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%Bz(i,j,k0+1) - this%Bz(i,j,k0)
    delta_plus  = this%Bz(i,j,k0+2) - this%Bz(i,j,k0+1)
    BzR = this%Bz(i,j,k0+1) - half * delta(delta_minus, delta_plus)
    ! energy
    delta_minus = this%e(i,j,k0) - this%e(i,j,k0-1)
    delta_plus  = this%e(i,j,k0+1) - this%e(i,j,k0)
    eL = this%e(i,j,k0) + half * delta(delta_minus, delta_plus)

    delta_minus = this%e(i,j,k0+1) - this%e(i,j,k0)
    delta_plus  = this%e(i,j,k0+2) - this%e(i,j,k0+1)
    eR = this%e(i,j,k0+1) - half * delta(delta_minus, delta_plus)

  end subroutine MINMOD_isolated
  !============================================================================!
  real(kind=8) function delta(delta_minus, delta_plus)
    implicit none
    real(kind=8), intent(in) :: delta_minus, delta_plus
    real(kind=8), parameter :: beta = 1.0d0

    if(delta_plus > zero) then
      delta = max(zero, min(beta * delta_minus, delta_plus), min(delta_minus, beta * delta_plus))
    else
      delta = min(zero, max(beta * delta_minus, delta_plus), max(delta_minus, beta * delta_plus))
    end if

  end function delta
  !============================================================================!
  subroutine boundary_isolated(this)
    implicit none
    class(hydro_base), intent(in out) :: this

    if(row_id.eq.0) this%u(:,1,:,:) = this%u(:,2,:,:)
    if(row_id.eq.row_nproc-1) this%u(:,domain%Nx,:,:) = this%u(:,domain%Nx-1,:,:)

    if(col_id.eq.0) this%u(:,:,1,:) = this%u(:,:,2,:)
    if(col_id.eq.col_nproc-1) this%u(:,:,domain%Ny,:) = this%u(:,:,domain%Ny-1,:)

    this%u(:,:,:,1) = this%u(:,:,:,2)
    this%u(:,:,:,domain%Nz) = this%u(:,:,:,domain%Nz-1)

  end subroutine boundary_isolated
  !============================================================================!

end module hydro_base_lib
