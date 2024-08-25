subroutine solve_Poisson
  use global_numbers
  use integrals
  use fft_lib, only: fft, ifft, sinfft, isinfft
  implicit none
  integer :: i, j, k
  real(kind=8) :: M, sponge_rc, sponge_delta, sponge_Vcoeff
  complex(kind=8), parameter :: i_unreal = (0.0d0, 1.0d0)
  complex(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: V_local

  gas%rho = max(gas%u(1,:,:,:), gas%floor)
  M = trapezium(gas%rho, domain%dx, domain%dy, domain%dz)

  if(gas%boundary.eq.'periodic') then

    V_local = fft(dcmplx(gas%rho - M / (domain%xmax - domain%xmin)**3))

    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=1, domain%Nz
          if(domain%k_square(i,j,k).ne.0.0d0)  then
            V_local(i,j,k) = - V_local(i,j,k) / domain%k_square(i,j,k)
          else
            V_local(i,j,k) = 0.0d0
          end if
        end do
      end do
    end do

    V_local = ifft(V_local);

  else if(gas%boundary.eq.'isolated') then
  
    V_local = sinfft(dcmplx(gas%rho))
  
    V_local(:,:,1) = 0.0d0
    do i=1, domain%Nx
      do j=1, domain%Ny
        do k=2, domain%Nz
          V_local(i,j,k) = 0.50d0 * domain%dz**2 * V_local(i,j,k) / (&
          dcos(pi * dble(i-1+coords(1)*domain%Nx) / dble(domain%Nz)) + &
          dcos(pi * dble(j-1+coords(2)*domain%Ny) / dble(domain%Nz)) + &
          dcos(pi * dble(k-1) / dble(domain%Nz)) - 3.0d0)
          if(coords(1).eq.0.and.i.eq.1) V_local(i,j,k) = 0.0d0
          if(coords(2).eq.0.and.j.eq.1) V_local(i,j,k) = 0.0d0
        end do
      end do
    end do
  
    V_local = isinfft(V_local)
  
  end if

  gas%V = dreal(V_local)

end subroutine solve_Poisson
!============================================================================!
