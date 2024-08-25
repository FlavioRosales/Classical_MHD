module global_numbers
  use mesh_lib
  use hydro_lib
  implicit none
  !
  ! numbers of points
  !
  integer :: Nx, Ny, Nz
  !
  ! domain
  !
  real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax, tmin, tmax
  !
  ! Courant number
  !
  real(kind=8) :: CFL
  !
  ! Save data every some steps
  !
  integer :: save_0d, save_1d, save_2d, save_3d
  !
  ! measured of cpu time
  !
  real(kind=8) :: t_star, t_cpu
  !
  ! Complete system
  !
  character(len=100) :: initial_data, file
  logical :: read_from_check_point = .false.
  character(len=20) :: boundary
  real(kind=8), allocatable, dimension(:) :: r
  !
  ! Euler system
  !
  type(hydro) :: gas
  character(len=20) :: EoS, flux_formula
  real(kind=8) :: floor, gamma
  real(kind=8) :: omega_star, R_star, M_star, K_poly
  real(kind=8) :: Mgas, Ugas, Kgas, Wgas, Egas, Qgas, rhomaxgas, pgas
  real(kind=8), dimension(2) :: pxgas, pygas, pzgas
  real(kind=8), dimension(3) :: xcmgas
  logical :: fourier_space_gas

end module global_numbers
