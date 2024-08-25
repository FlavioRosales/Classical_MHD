subroutine initial_conditions
  use global_numbers
  implicit none
  
  if(initial_data.eq.'equilibrium') then
    call gas%equilibrium(xcm = [0.0d0, 0.0d0, 0.0d0], vcm = [0.0d0, 0.0d0, 0.0d0], M_star = M_star, R_star = R_star)
  else if(initial_data.eq.'binary') then
    call gas%binary(xcm_1 = [-10.0d0, +10.0d0, 0.0d0], vcm_1 = [0.050d0, 0.0d0, 0.0d0], &
    M_star_1 = M_star, R_star_1 = R_star, R_star_2 = R_star/2.0d0, MR = 1.0d0)
  end if
  call solve_Poisson

end subroutine initial_conditions
!==============================================================================!
! subroutine dark_matter
!   use global_numbers
!   use read_data_lib
!   use interpolation_lib
!   use integrals
!   implicit none
!   real(kind=8) :: sigmad0, rhob0, lambda, M200, M210, s
!   real(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: rho200, rho210

!   rho200 = abs(wave(filename = filename_haloformed, l=0, m=0, lambda = 1.0d0))**2
!   rho210 = abs(wave(filename = filename_haloformed, l=1, m=0, lambda = 1.0d0))**2
!   M200 = trapezium(rho200, domain%dx, domain%dy, domain%dz)
!   M210 = trapezium(rho210, domain%dx, domain%dy, domain%dz)
!   if(rank.eq.master) print*, 'M200 / M210 = ', M200 / M210
!   rho =  rho200 + rho210

! contains

!   function wave(filename, l, m, lambda) result(psi)
!     implicit none
!     !
!     ! input
!     !
!     character(len=*), intent(in) :: filename
!     integer, intent(in) :: l, m
!     real(kind=8), intent(in) :: lambda
!     !
!     ! output
!     !
!     complex(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: psi
!     !
!     ! internal
!     !
!     integer :: i, j, k, n
!     type(read_data) :: data
!     real(kind=8) :: atom
!     complex(kind=8) :: i_unreal = (0.0d0, 1.0d0)
!     real(kind=8), allocatable, dimension(:) :: r, psi_r
!     complex(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: Ylm

!     data = read_data(columns_number = 14, filename = '../../../IO/input/multistates/' // filename)

!     allocate(r(0:data%lines_number), psi_r(0:data%lines_number))

!     r = data%data(1,:) / lambda

!     if(l.eq.0) then

!       Ylm = 1.0d0
!       atom = 1.0d0

!       psi_r = data%data(2,:) * lambda**2

!     end if

!     if(l.eq.1) then

!       if(m.eq.-1) then

!         Ylm = + dsqrt(3.0d0/8.0d0/pi) * (domain%x - i_unreal * domain%y)

!       else if(m.eq.0) then

!         atom = 3.544907701811030d0
!         Ylm = dsqrt(3.0D0/4.0D0/pi) * domain%z

!       else if(m.eq.1) then

!         Ylm = - dsqrt(3.0d0/8.0d0/pi) * (domain%x + i_unreal * domain%y)

!       end if

!       psi_r = data%data(4,:) * lambda**3

!     end if


!     do i=1, domain%Nx
!       do j=1, domain%Ny
!         do k=1, domain%Nz

!           if(domain%r(i,j,k) < r(data%lines_number)) then
!             n = min(max(int((domain%r(i,j,k) - r(0)) / (r(1) - r(0))), 1), data%lines_number-1)
!             psi(i,j,k) = Lagrange(x = [r(n-1), r(n), r(n+1)], f = [psi_r(n-1), psi_r(n), psi_r(n+1)], point = domain%r(i,j,k))
!           else
!             psi(i,j,k) = 0.0d0
!           end if

!         end do
!       end do
!     end do

!     psi = atom * psi * Ylm !/ dsqrt(4.0d0 * pi)

!   end function wave

! end subroutine dark_matter
! !==============================================================================!
