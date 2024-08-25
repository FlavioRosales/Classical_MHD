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

  end type hydro_initial_data

contains

  !=============================================================================!
  !=============================================================================!
  

end module hydro_initial_data_lib
