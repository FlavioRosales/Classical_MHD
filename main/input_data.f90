subroutine input_data
  use strings_lib
  use global_numbers
  implicit none
  character(len=40) :: input_file
  type(string) :: filename
  type(string) :: filename_split(2)
  type(string) :: output_folder

  !
  ! create the list of parameters
  !
  namelist/input/ Nx, Ny, Nz, &
  xmin, xmax, ymin, ymax, zmin, zmax, tmin, tmax, CFL, &
  save_0d, save_1d, save_2d, save_3d, &
  omega_star, R_star, M_star, gamma, floor, K_poly, &
  EoS, flux_formula, initial_data, file, boundary, read_from_check_point
  !
  ! Default Values
  !
  Nx = 128
  Ny = 128
  Nz = 128
  xmin = -50.0d0
  xmax = +50.0d0
  ymin = -50.0d0
  ymax = +50.0d0
  zmin = -50.0d0
  zmax = +50.0d0
  tmin = 0.0d0
  tmax = 100.0d0
  CFL = 0.250d0
  save_0d = 1
  save_1d = 10
  save_2d = 100
  save_3d = 100
  omega_star = 0.0d0
  R_star = 5.0d0
  M_star = 10.0d0
  gamma = 5.0d0 / 3.0d0
  floor = 1.0e-10
  K_poly = 10.0d0
  EoS = 'ideal gas'
  flux_formula = 'HLLE'
  initial_data = 'equilibrium'
  boundary = 'periodic'
  read_from_check_point = .false.
  file = 'data_IC_hidro_z100'
  !
  ! read the name of the parameter file
  !
  boundary = trim(boundary)
  call get_command_argument(1, input_file); input_file = trim(input_file)

  if(rank.eq.master) print*, ' '
  if(rank.eq.master) print*, 'Reading input parameters from file: ' // input_file
  !
  ! now the list of parameters is read
  !
  open(100, file = input_file, action = 'read')
    read(100, nml = input)
  close(100)
  !
  ! eliminate the spaces
  !
  filename = remove(input_file, '')
  !
  ! remove the extension ".something"
  !
  filename_split = filename%split('.')
  output_folder = filename_split(1) + '/'
  !
  ! create the folder where the output data is saved.
  !
  call system('mkdir -p '//trim(output_folder%string_data))
  !
  ! change the path to the folder where the data is saved
  !
  call CHDIR(output_folder%string_data)

end subroutine input_data
