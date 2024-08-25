subroutine initialize_variables
  use global_numbers
  use hdf5_lib
  use strings_lib
  implicit none

  integer :: i
  character(len=100) :: F; F = '(A, F10.3, A, F10.3, A, ES14.7, A)'

  if(rank.eq.master) print*, ' '
  if(rank.eq.master) print*, 'Initializing system variables...'

  !
  ! Initialize the mes variable
  !
  domain = mesh(&
  xmin = xmin, xmax = xmax, Nx = Nx, &
  ymin = ymin, ymax = ymax, Ny = Ny, &
  zmin = zmin, zmax = zmax, Nz = Nz, &
  tmin = tmin, tmax = tmax, CFL = CFL)

  call domain%save_indexes(x_save = 0.0d0, y_save = 0.0d0, z_save = 0.0d0)
  !
  ! Initialize hydro variables
  !
  gas  = hydro(K_poly = K_poly, gamma = gamma, floor = floor, EoS = EoS, flux_formula = flux_formula, boundary = trim(boundary))
  !
  ! hdf5 files: gas and bec
  !
  !
  ! create files
  !
  call create_hdf5_file('rho.h5' )
  call create_hdf5_file('v.h5' )
  !
  ! create groups
  !
  call create_hdf5_group('rho.h5' , '/refinement_1')
  call create_hdf5_group('v.h5' , '/refinement_1')
  !
  ! save grid
  !
  call write_hdf5('rho.h5' , '/refinement_1', 'Xcoord', domain%x)
  call write_hdf5('rho.h5' , '/refinement_1', 'Ycoord', domain%y)
  call write_hdf5('rho.h5' , '/refinement_1', 'Zcoord', domain%z)
  call write_hdf5('v.h5' , '/refinement_1', 'Xcoord', domain%x)
  call write_hdf5('v.h5' , '/refinement_1', 'Ycoord', domain%y)
  call write_hdf5('v.h5' , '/refinement_1', 'Zcoord', domain%z)

  !
  ! Numerical domain:
  !
  if(rank.eq.master) then
    print*, ' '
    write(*,*) 'Numerical domain: '
    print*, ' '
    write(*,*) '|===========================================================|'
    write(*,F) ' | xmin = ', xmin, '| xmax = ', xmax, '| dx = ', domain%dx, ' |'
    write(*,F) ' | ymin = ', ymin, '| ymax = ', ymax, '| dy = ', domain%dy, ' |'
    write(*,F) ' | zmin = ', zmin, '| zmax = ', zmax, '| dz = ', domain%dz, ' |'
    write(*,F) ' | tmin = ', tmin, '| tmax = ', tmax, '| dt = ', domain%dt, ' |'
    write(*,*) '|===========================================================|'
  end if

end subroutine initialize_variables
