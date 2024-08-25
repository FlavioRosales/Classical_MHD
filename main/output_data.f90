subroutine output_data(index)
  use save_ASCII
  use global_numbers
  use strings_lib
  use hdf5_lib
  implicit none
  integer, intent(in) :: index
  real(kind=8), dimension(3) :: xcm
  integer :: i
  character(len=31) :: F; F = '(A, I10, A, F10.3, A, F10.3, A)'

  !
  ! Save data along of time axis
  !
  if(mod(index, save_0d).eq.0) then

    call overall_diagnosis

    call save_t(domain%t, rhomaxgas, 'rhomax')
    call save_t(domain%t, mgas, 'M')
    call save_t(domain%t, Kgas, 'K')
    call save_t(domain%t, Wgas, 'W')
    call save_t(domain%t, Ugas, 'U')
    call save_t(domain%t, Egas, 'E')
    call save_t(domain%t, Qgas, 'Q')
    call save_t(domain%t, pxgas, 'px')
    call save_t(domain%t, pygas, 'py')
    call save_t(domain%t, pzgas, 'pz')
    call save_t(domain%t,  pgas, 'p')
    call save_t(domain%t, xcmgas, 'xcm')

  end if
  !
  ! Save data along of the x-axis
  !
  if(domain%save_yz .and. mod(index, save_1d).eq.0) then

    call save_x(domain%x(:, domain%j_save, domain%k_save), &
    gas%rho(:, domain%j_save, domain%k_save), 'rho')
    call save_x(domain%x(:, domain%j_save, domain%k_save), &
    gas%e(:, domain%j_save, domain%k_save), 'e')
    call save_x(domain%x(:, domain%j_save, domain%k_save), &
    gas%vx(:, domain%j_save, domain%k_save), 'vx')
    call save_x(domain%x(:, domain%j_save, domain%k_save), &
    gas%V(:, domain%j_save, domain%k_save), 'V')
    call save_x(domain%x(:, domain%j_save, domain%k_save), &
    gas%p(:, domain%j_save, domain%k_save) / gas%rho(:, domain%j_save, domain%k_save)**gas%gamma , 'K')

  end if
  !
  ! Save data along of the y-axis
  !
  if(domain%save_xz .and. mod(index, save_1d).eq.0) then

    call save_y(domain%y(domain%i_save, :, domain%k_save), &
    gas%rho(domain%i_save, :, domain%k_save), 'rho')
    call save_y(domain%y(domain%i_save, :, domain%k_save), &
    gas%e(domain%i_save, :, domain%k_save), 'e')
    call save_y(domain%y(domain%i_save, :, domain%k_save), &
    gas%vy(domain%i_save, :, domain%k_save), 'vy')

  end if
  !
  ! Save data along of the z-axis
  !
  if(domain%save_xy .and. mod(index, save_1d).eq.0) then

    call save_z(domain%z(domain%i_save, domain%j_save, :), &
    gas%rho(domain%i_save, domain%j_save, :), 'rho')
    call save_z(domain%z(domain%i_save, domain%j_save, :), &
    gas%e(domain%i_save, domain%j_save, :), 'e')
    call save_z(domain%z(domain%i_save, domain%j_save, :), &
    gas%vz(domain%i_save, domain%j_save, :), 'vz')

  end if
  !
  ! Save data at the plane xy
  !
  if(domain%save_z .and. mod(index, save_2d).eq.0) then

    call save_xy(domain%x(:,:,domain%k_save), domain%y(:,:,domain%k_save), &
    gas%rho(:,:,domain%k_save), 'rho')
    call save_xy(domain%x(:,:,domain%k_save), domain%y(:,:,domain%k_save), &
    gas%e(:,:,domain%k_save), 'e')
    call save_xy_vect(domain%x(:,:,domain%k_save), domain%y(:,:,domain%k_save), &
    gas%vx(:,:,domain%k_save), gas%vy(:,:,domain%k_save), 'vxvy')
    call save_xy(domain%x(:,:,domain%k_save), domain%y(:,:,domain%k_save), &
    gas%p(:,:,domain%k_save) / gas%rho(:,:,domain%k_save)**gas%gamma, 'K')

  end if
  !
  ! Save data at the plane xz
  !
  if(domain%save_y .and. mod(index, save_2d).eq.0) then

    call save_xz(domain%x(:,domain%j_save,:), domain%z(:,domain%j_save,:), &
    gas%rho(:,domain%j_save,:), 'rho')
    call save_xz(domain%x(:,domain%j_save,:), domain%z(:,domain%j_save,:), &
    gas%e(:,domain%j_save,:), 'e')
    call save_xz_vect(domain%x(:,domain%j_save,:), domain%z(:,domain%j_save,:), &
    gas%vx(:,domain%j_save,:), gas%vz(:,domain%j_save,:), 'vxvz')

  end if
  !
  ! Save data in 3d space
  !
  if(mod(index, save_3d).eq.0) then

    call write_hdf5('rho.h5' , '/refinement_1', 'rho_' // char_integer(index/save_3d), gas%rho)
    call write_hdf5('v.h5' , '/refinement_1', 'vx_' // char_integer(index/save_3d), gas%vx)
    call write_hdf5('v.h5' , '/refinement_1', 'vy_' // char_integer(index/save_3d), gas%vy)
    call write_hdf5('v.h5' , '/refinement_1', 'vz_' // char_integer(index/save_3d), gas%vz)
    call gas%write_check_point(name =  'gas', time = domain%t)

  end if


  if(rank.eq.master .and. mod(index, save_2d).eq.0) then
    write(*,*) '|======================================|'
    write(*,*) '|  index     |   time     |  cpu time  |'
  end if

  if(rank.eq.master .and. mod(index,save_1d).eq.0) then
    write(*,*) '|======================================|'
    write(*,F) ' | ', index, ' | ', domain%t, ' | ', t_cpu, ' |'
  end if



end subroutine output_data
