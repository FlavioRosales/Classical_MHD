module hydro_initial_data_lib
  use hydro_diagnostics_lib
  implicit none

  integer, private :: i_star
  integer, parameter, private :: Nxi = 10000
  real(kind=8), dimension(0:Nxi), private :: xi, r, rho
  real(kind=8), dimension(2,0:Nxi), private :: u
  real(kind=8), parameter, private :: G = 1.0d0 / (4.0d0 * pi)
  real(kind=8), private :: n, xi_star, number_N, alpha, rhoc, rhoav, psi_star, K_poly

  type, extends(hydro_diagnostics) :: hydro_initial_data
  contains

    procedure :: boost
    procedure :: equilibrium
    procedure :: binary
    procedure :: multiple_clouds
    procedure :: equilibrium_sphwavegas

  end type hydro_initial_data

contains

  !=============================================================================!
  subroutine spherical_solution(M_star, R_star)
    use ODE
    implicit none
    real(kind=8) :: M_star, R_star

    xi = linspace(min = 0.0d0, max = 100.0d0, npoints = Nxi)
    xi_star = 0.0d0
    u = solve(initial_conditions = [1.0d0, 0.0d0], rhs = Lane_Emden, domain = xi)
    u(1,:) = max(u(1,:), zero)
    i_star = int(xi_star / (xi(1) - xi(0)))
    psi_star = abs(u(2,i_star))
    number_N = (4.0d0 * pi)**(1.0d0 / n) / (n + 1.0d0) * psi_star**((1.0d0 - n) / n) * xi_star**((n - 3.0d0) / n)
    rhoav = 3.0d0 * M_star / (4.0d0 * pi * R_star**3)
    rhoc = rhoav * xi_star**3 / (3.0d0 * psi_star)
    K_poly = G * number_N * R_star**((3.0d0 - n)/n) * M_star**((n - 1.0d0)/n)
    alpha = sqrt(K_poly * (n + one) * rhoc**((one - n)/n) / (4.0d0 * pi * G))
    r = alpha * xi
    rho = rhoc * u(1,:)**n

  end subroutine spherical_solution
  !=============================================================================!
  function Lane_Emden(xi_local, u_local) result(rhs)
    implicit none
    real(kind=8), intent(in) :: xi_local, u_local(:)
    real(kind=8) :: rhs(size(u_local))

    real(kind=8) :: theta, psi

    theta = u_local(1)
    psi   = u_local(2) ! where psi if defined by \psi = \xi**2  \theta_{\xi}

    if(theta > 0.0d0) then

      if(xi_local.eq.0.0d0) then
        rhs = 0.0d0
      else
        rhs(1) = psi / xi_local**2
        rhs(2) = - xi_local**2 * theta**n
      end if

    else

      theta = 0.0d0
      if(xi_star.eq.0.0d0) xi_star = xi_local
      rhs = 0.0d0

    end if

  end function Lane_Emden
  !=============================================================================!
  function interpolate(r, ur, xcm) result(u)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(:), intent(in) :: r, ur, xcm
    !
    ! output
    !
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: u
    !
    ! internal
    !
    integer :: i, j, k, l, i1
    real(kind=8) :: ro

    do i = 1, domain%Nx
      do j = 1, domain%Ny
        do k = 1, domain%Nz

          ro = sqrt((domain%x(i,j,k) - xcm(1))**2 + (domain%y(i,j,k) - xcm(2))**2 + (domain%z(i,j,k) - xcm(3))**2)

          i1 = int((ro - r(1))/(r(2) - r(1))) + 1

          if(i1 < size(r)) then
            u(i,j,k) = ur(i1) + (ur(i1+1) - ur(i1)) / (r(i1+1) - r(i1)) * (ro - r(i1))
          else
            u(i,j,k) = zero
          end if

        end do
      end do
    end do

  end function interpolate
  !=============================================================================!
  function boost(this, v0) result(v)
    implicit none
    !
    ! input
    !
    class(hydro_initial_data), intent(in out) :: this
    real(kind=8), intent(in) :: v0
    !
    ! output
    !
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: v
    !
    ! internal
    !
    integer ::i, j, k, l

    do i = 1, domain%Nx
      do j = 1, domain%Ny
        do k = 1, domain%Nz

          if(this%rho(i,j,k) > this%floor) then
            v(i,j,k) = v0
          else
            v(i,j,k) = zero
          end if

        end do
      end do
    end do

  end function boost
  !=============================================================================!
  subroutine equilibrium(this, xcm, vcm, M_star, R_star)
    implicit none
    !
    ! input
    !
    class(hydro_initial_data), intent(in out) :: this
    real(kind=8), intent(in) :: M_star, R_star
    real(kind=8), dimension(3), intent(in) :: xcm, vcm

    integer :: i, j, k, l
    real(kind=8) :: r0

    n = one / (this%gamma - one)
    call spherical_solution(M_star = M_star, R_star = R_star)
    this%K_poly = K_poly

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          r0 = dsqrt((domain%x(i,j,k) - xcm(1))**2 + (domain%y(i,j,k) - xcm(2))**2 + (domain%z(i,j,k) - xcm(3))**2)
          l = int((r0 - r(0)) / (r(1) - r(0)))

          if(l<Nxi) then
            this%rho(i,j,k) = rho(l) + (rho(l+1) - rho(l)) * (r0 - r(l)) / (r(l+1) - r(l))
          else
            this%rho(i,j,k) = zero
          end if

          this%rho(i,j,k) = max(this%rho(i,j,k), this%floor)

          if(this%rho(i,j,k)>this%floor) then
            this%vx(i,j,k) = vcm(1)
            this%vy(i,j,k) = vcm(2)
            this%vz(i,j,k) = vcm(3)
          else
            this%vx(i,j,k) = zero
            this%vy(i,j,k) = zero
            this%vz(i,j,k) = zero
          end if

        end do
      end do
    end do

    if(this%EoS.eq.'ideal gas' .or. this%EoS.eq.'poly') then
      this%p = this%K_poly * this%rho**this%gamma
      this%e = this%p / this%rho / (this%gamma - one)
    else if(this%EoS.eq.'dust') then
      this%p = zero
      this%e = zero
    end if

    this%u(1,:,:,:) = this%rho
    this%u(2,:,:,:) = this%rho * this%vx
    this%u(3,:,:,:) = this%rho * this%vy
    this%u(4,:,:,:) = this%rho * this%vz
    this%u(5,:,:,:) = this%rho * (this%e + half * (this%vx**2 + this%vy**2 + this%vz**2))

  end subroutine equilibrium
  !=============================================================================!
  subroutine binary(this, xcm_1, vcm_1, M_star_1, R_star_1, R_star_2, MR)
    implicit none
    !
    ! input
    !
    class(hydro_initial_data), intent(in out) :: this
    real(kind=8), intent(in) :: M_star_1, R_star_1, R_star_2, MR
    real(kind=8), dimension(3), intent(in) :: xcm_1, vcm_1

    integer :: i, j, k, l
    real(kind=8) :: M_star_2, K_poly_1, K_poly_2, r0
    real(kind=8), dimension(3) :: xcm_2, vcm_2
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: rho1, rho2, p1, p2


    M_star_2 = M_star_1 / MR
    xcm_2 = - MR * xcm_1
    vcm_2 = - MR * vcm_1

    n = one / (this%gamma - one)
    call spherical_solution(M_star = M_star_1, R_star = R_star_1)
    K_poly_1 = K_poly

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          r0 = dsqrt((domain%x(i,j,k) - xcm_1(1))**2 + (domain%y(i,j,k) - xcm_1(2))**2 + (domain%z(i,j,k) - xcm_1(3))**2)
          l = int((r0 - r(0)) / (r(1) - r(0)))

          if(l<Nxi) then
            rho1(i,j,k) = rho(l) + (rho(l+1) - rho(l)) * (r0 - r(l)) / (r(l+1) - r(l))
          else
            rho1(i,j,k) = zero
          end if

          rho1(i,j,k) = max(rho1(i,j,k), this%floor)

        end do
      end do
    end do

    call spherical_solution(M_star = M_star_2, R_star = R_star_2)
    K_poly_2 = K_poly

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          r0 = dsqrt((domain%x(i,j,k) - xcm_2(1))**2 + (domain%y(i,j,k) - xcm_2(2))**2 + (domain%z(i,j,k) - xcm_2(3))**2)
          l = int((r0 - r(0)) / (r(1) - r(0)))

          if(l<Nxi) then
            rho2(i,j,k) = rho(l) + (rho(l+1) - rho(l)) * (r0 - r(l)) / (r(l+1) - r(l))
          else
            rho2(i,j,k) = zero
          end if

          rho2(i,j,k) = max(rho2(i,j,k), this%floor)

        end do
      end do
    end do

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz

          if(rho1(i,j,k)>this%floor) then
            this%vx(i,j,k) = vcm_1(1)
            this%vy(i,j,k) = vcm_1(2)
            this%vz(i,j,k) = vcm_1(3)
          end if

          if(rho2(i,j,k)>this%floor) then
            this%vx(i,j,k) = vcm_2(1)
            this%vy(i,j,k) = vcm_2(2)
            this%vz(i,j,k) = vcm_2(3)
          end if

        end do
      end do
    end do

    p1 = K_poly_1 * rho1**this%gamma
    p2 = K_poly_2 * rho2**this%gamma

    if(this%EoS.eq.'ideal gas' .or. this%EoS.eq.'poly') then
      this%rho = max(rho1 + rho2, this%floor)
      this%p = p1 + p2
      this%e = this%p / this%rho / (this%gamma - one)
    else if(this%EoS.eq.'dust') then
      this%rho = max(rho1 + rho2, this%floor)
      this%p = zero
      this%e = zero
    end if

    this%u(1,:,:,:) = this%rho
    this%u(2,:,:,:) = this%rho * this%vx
    this%u(3,:,:,:) = this%rho * this%vy
    this%u(4,:,:,:) = this%rho * this%vz
    this%u(5,:,:,:) = this%rho * (this%e + half * (this%vx**2 + this%vy**2 + this%vz**2))

  end subroutine binary
  !=============================================================================!
  subroutine multiple_clouds(this, n_cores, M_total, fixed_seed)
    implicit none
    !
    ! input
    !
    class(hydro_initial_data), intent(in out) :: this
    real(kind=8), intent(in) :: M_total
    integer, intent(in) :: n_cores
    logical, intent(in) :: fixed_seed
    !
    ! internal
    !
    integer :: i, j, k, l, n
    real(kind=8) :: min, max
    real(kind=8), dimension(n_cores) :: &
    x_local, y_local, z_local, vx_local, vy_local, vz_local, mass_local, r_local
    real(kind=8), dimension(n_cores) :: x, y, z, vx, vy, vz, mass, lambda, r
    real(kind=8) :: M
    real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: rho3d, vr, r3d

    !
    ! Initialise random number generator.
    !
    if(.not.fixed_seed) then
      call random_seed(size=n)              ! Get size of seed array.
      call random_seed(put=urandom_seed(n)) ! Put seed array into PRNG.
    end if
    !
    ! Generate random numbers
    !
    call random_number(x_local)
    call random_number(y_local)
    call random_number(z_local)
    call random_number(mass_local)
    call random_number(r_local)

    call MPI_ALLREDUCE(x_local, x, n_cores, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(y_local, y, n_cores, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(z_local, z, n_cores, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(mass_local, mass, n_cores, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(r_local, r, n_cores, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    x = x - minval(x); x = x / maxval(x)
    y = y - minval(y); y = y / maxval(y)
    z = z - minval(z); z = z / maxval(z)
    mass = mass - minval(mass); mass = mass / maxval(mass)
    r = r - minval(r); r = r / maxval(r)

    x = (domain%xmax - domain%xmin) * x + domain%xmin
    y = (domain%ymax - domain%ymin) * y + domain%ymin
    z = (domain%zmax - domain%zmin) * z + domain%zmin

    mass = mass * 0.80d0 + 0.20d0

    r = r * 10.0d0 + 2.0d0

    M = 0.0d0
    do i=1, n_cores
      M = M + mass(i)
    end do
    mass = mass * M_total / M

     x(n_cores) = 0.0d0;  y(n_cores) = 0.0d0;  z(n_cores) = 0.0d0

    do i=1, n_cores-1
      x(n_cores) = x(n_cores) - mass(i) * x(i)
      y(n_cores) = y(n_cores) - mass(i) * y(i)
      z(n_cores) = z(n_cores) - mass(i) * z(i)
    end do

    x(n_cores) = x(n_cores) / mass(n_cores)
    y(n_cores) = y(n_cores) / mass(n_cores)
    z(n_cores) = z(n_cores) / mass(n_cores)

    x = 0.30d0 * (domain%xmax - domain%xmin) * x / maxval(abs(x))
    y = 0.30d0 * (domain%ymax - domain%ymin) * y / maxval(abs(y))
    z = 0.30d0 * (domain%zmax - domain%zmin) * z / maxval(abs(z))

    n = one / (this%gamma - one)

    this%rho = zero
    this%  p = zero
    this% vx = zero
    this% vy = zero
    this% vz = zero

    do i=1, n_cores

      if(this%EoS.eq.'ideal gas' .or. this%EoS.eq.'poly') then
        call spherical_solution(M_star = mass(i), R_star = r(i))
        this%K_poly = K_poly

        rho3d = max(interpolate(r = r, ur = rho, xcm = [x(i), y(i), z(i)]), zero)
        this%rho = this%rho + rho3d
        this%  p = this%p + this%K_poly * max(rho3d, this%floor)

      else if(this%EoS.eq.'dust') then


        r3d = sqrt((domain%x - x(i))**2 + (domain%y - y(i))**2 + (domain%z - z(i))**2)
        this%rho = this%rho + 3.0d0 * mass(i) / (4.0d0 * pi * r(i)**2) * (1.0d0 + r3d**2/r(i)**2)**(-2.50d0)
        vr = sqrt(2.0d0 * G * mass(i) / (6.0d0 * sqrt(r3d**2 + r(i)**2)))

        do l=1, domain%Nx
          do j=1, domain%Ny
            do k=1, domain%Nz

                if(r3d(l,j,k).ne.0.0d0) then
                  this%vx(l,j,k) = this%vx(l,j,k) + vr(l,j,k) * (domain%x(l,j,k) - x(i)) / r3d(l,j,k)
                  this%vy(l,j,k) = this%vy(l,j,k) + vr(l,j,k) * (domain%y(l,j,k) - y(i)) / r3d(l,j,k)
                  this%vz(l,j,k) = this%vz(l,j,k) + vr(l,j,k) * (domain%z(l,j,k) - z(i)) / r3d(l,j,k)
                end if

            end do
          end do
        end do


      end if

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    end do

    this%rho = max(this%rho, this%floor)
    this%e  = this%p / this%rho / (this%gamma - 1.0d0)

    this%u(1,:,:,:) = this%rho
    this%u(2,:,:,:) = this%rho * this%vx
    this%u(3,:,:,:) = this%rho * this%vy
    this%u(4,:,:,:) = this%rho * this%vz
    this%u(5,:,:,:) = this%rho * (this%e + half * (this%vx**2 + this%vy**2 + this%vz**2))

  end subroutine multiple_clouds
  !=============================================================================!
  function urandom_seed(n, stat) result(seed)
      !! Returns a seed array filled with random values from `/dev/urandom`.
      integer, intent(in)            :: n
      integer, intent(out), optional :: stat
      integer                        :: seed(n)
      integer                        :: fu, rc

      open (access='stream', action='read', file='/dev/urandom', &
            form='unformatted', iostat=rc, newunit=fu)
      if (present(stat)) stat = rc
      if (rc == 0) read (fu) seed
      close (fu)

  end function urandom_seed
  !============================================================================!
  subroutine equilibrium_sphwavegas(this, name)
    use read_data_lib
    use interpolation_lib
    implicit none
    class(hydro_initial_data), intent(in out) :: this
    character(len=*), intent(in) :: name
    !
    ! internal
    !
    integer :: i, j, k, l
    type(read_data) :: data
    real(kind=8) :: r

    this%K_poly = 10.0d0

    data = read_data(filename = '../../../IO/input/' // name, columns_number = 3)

    if(rank.eq.master) print*, '     Interpolating a BEC core...'

    do i = 1, domain%Nx
      do j = 1, domain%Ny
        do k = 1, domain%Nz

          r = sqrt(domain%x(i,j,k)**2 + domain%y(i,j,k)**2 + domain%z(i,j,k)**2)

          l = int((r - data%data(1,0)) / (data%data(1,1) - data%data(1,0)))

          if(l < data%lines_number) then
            this%rho(i,j,k) =max(Lagrange(point = r, &
            x = [data%data(1,l), data%data(1,l+1)], &
            f = [data%data(3,l), data%data(3,l+1)]), this%floor)
          else
            this%rho(i,j,k) = this%floor
          end if

        end do
      end do
    end do

    this%p = this%K_poly * this%rho**this%gamma
    this%e = this%p / this%rho / (this%gamma - one)

    this%u(1,:,:,:) = this%rho
    this%u(2,:,:,:) = zero
    this%u(3,:,:,:) = zero
    this%u(4,:,:,:) = zero
    this%u(5,:,:,:) = this%rho * this%e

  end subroutine equilibrium_sphwavegas
  !==============================================================================!

end module hydro_initial_data_lib
