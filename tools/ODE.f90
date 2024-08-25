module ODE
  implicit none

contains
  !=============================================================================!
  function linspace(min, max, npoints) result(domain)
    implicit none
    integer, intent(in) :: npoints
    real(kind=8), intent(in) :: min, max
    real(kind=8), dimension(0:npoints) :: domain
    integer :: i
    real(kind=8) :: h

    h = (max - min) / dble(npoints)

    do i=0, npoints
      domain(i) = min + dble(i) * h
    end do

  end function linspace
  !=============================================================================!
  function solve(initial_conditions, domain, rhs) result(solution)
    implicit none
    real(kind=8), dimension(:) :: initial_conditions, domain
    interface
      function rhs(t, x) result(dotx)
        real(kind=8), intent(in) :: t, x(:)
        real(kind=8), dimension(size(x)) :: dotx
      end function rhs
    end interface
    real(kind=8), dimension(size(initial_conditions), size(domain)) :: solution
    integer :: j
    real(kind=8) :: h
    real(kind=8), dimension(size(initial_conditions)) :: k1, k2, k3, k4

    solution(:,1) = initial_conditions

    do j=2, size(domain)

      h = domain(j) - domain(j-1)

      k1 = rhs(domain(j-1), solution(:,j-1)) * h
      k2 = rhs(domain(j-1) + 0.50d0 * h, solution(:,j-1) + 0.50d0 * k1) * h
      k3 = rhs(domain(j-1) + 0.50d0 * h, solution(:,j-1) + 0.50d0 * k2) * h
      k4 = rhs(domain(j-1) + h, solution(:,j-1) + k3) * h

      solution(:,j) = solution(:,j-1) + (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4) / 6.0d0

    end do

  end function solve
  !=============================================================================!
  real(kind=8) function linear_interpolate(f, domain, point)
    implicit none
    real(kind=8), intent(in) :: f(:), domain(:), point
    integer i

    do i=1, size(f)-1

      if(domain(i).le.point .and. point.le.domain(i+1)) then
        linear_interpolate = f(i) + (point - domain(i)) * (f(i+1) - f(i)) / (domain(i+1) - domain(i))
        return
      end if

    end do

  end function linear_interpolate
  !=============================================================================!

end module ODE
