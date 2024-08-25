subroutine system_evolve
  use global_numbers
  implicit none

  integer :: j, rk
  real(kind=8), dimension(5, domain%Nx, domain%Ny, domain%Nz) :: k1, k2, k3, k4

  if(rank.eq.master) print*, ''
  if(rank.eq.master) print*, 'Creating the initial conditions...'

  t_star = MPI_Wtime()

  call initial_conditions

  if(rank.eq.master) print*, ''
  if(rank.eq.master) print*, 'Beginning of evolution...'


  call output_data(0)

  t_cpu = MPI_Wtime() - t_star

  do j=1, domain%Nt

    gas%u_p   = gas%u
    domain%t_p = domain%t

    do rk=1, 4

      call solve_Poisson
      !
      ! Evolve: GAS
      !
      if(rk.eq.1) then
        k1 = gas%rhs_hydro() * domain%dt
        gas%u = gas%u_p + 0.50d0 * k1
        domain%t = domain%t_p + 0.50d0 * domain%dt
      else if(rk.eq.2) then
        k2 = gas%rhs_hydro() * domain%dt
        gas%u = gas%u_p + 0.50d0 * k2
        domain%t = domain%t_p + 0.50d0 * domain%dt
      else if(rk.eq.3) then
        k3 = gas%rhs_hydro() * domain%dt
        gas%u = gas%u_p + k3
        domain%t = domain%t_p + domain%dt
      else
        k4 = gas%rhs_hydro() * domain%dt
        gas%u = gas%u_p + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0
        domain%t = domain%t_p + domain%dt
      end if

      if(gas%boundary.eq.'isolated') call gas%boundary_isolated

    end do

    call output_data(j)
    t_cpu = MPI_Wtime() - t_star

  end do

end subroutine system_evolve
