module interpolation_lib
  implicit none

  interface Lagrange
    module procedure :: Lagrange_1d, Lagrange_2d, Lagrange_3d, Lagrange_3d_array
  end interface Lagrange

contains
  !============================================================================!
  ! Lagrange interpolation in 1D.
  !============================================================================!
  function Lagrange_1d(x, f, point) result(solution)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: point
    real(kind=8), dimension(:), intent(in) :: x, f
    !
    ! output
    !
    real(kind=8) :: solution
    !
    ! internal
    !
    integer :: i, j
    real(kind=8), dimension(size(x)) :: l

    solution = 0.0d0
    l  = 1.0d0

    do j=1, size(x)

      if(point.ne.x(j)) then

        do i=1, size(x)
          if(i.ne.j) l(j) = l(j) * (point - x(i)) / (x(j) - x(i))
        end do

        solution = solution + f(j) * l(j)

      else

        solution = f(j)

      end if

    end do

  end function Lagrange_1d
  !============================================================================!
  ! Lagrange interpolation in 2D.
  !============================================================================!
  function Lagrange_2d(x, y, f, point) result(solution)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(2), intent(in) :: point
    real(kind=8), dimension(:,:), intent(in) :: x, y, f
    !
    ! output
    !
    real(kind=8) :: solution
    !
    ! internal
    !
    integer :: i
    real(kind=8), dimension(size(y,2)) :: fy

    !
    ! we first interpolate along the x-axis.
    !
    do i=1, size(y,2)
      fy(i) = Lagrange_1d(x = x(:,1), f = f(:,i), point = point(1))
    end do
    !
    ! Now we interpolate along the y-axis.
    !
    solution = Lagrange_1d(x = y(1,:), f = fy, point = point(2))

  end function Lagrange_2d
  !============================================================================!
  ! Lagrange interpolation in 3D.
  !============================================================================!
  function Lagrange_3d(x, y, z, f, point) result(solution)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(3), intent(in) :: point
    real(kind=8), dimension(:,:,:), intent(in) :: x, y, z, f
    !
    ! output
    !
    real(kind=8) :: solution
    !
    ! internal
    !
    integer :: i
    real(kind=8), dimension(size(z,3)) :: fz

    !
    ! we first interpolate along the x-axis.
    !
    do i=1, size(z,3)
      fz(i) = Lagrange_2d(x = x(:,:,1), y = y(:,:,1), f = f(:,:,i), point = point(1:2))
    end do
    !
    ! Now we interpolate along the y-axis.
    !
    solution = Lagrange_1d(x = z(1,1,:), f = fz, point = point(3))

  end function Lagrange_3d
  !============================================================================!
  ! Lagrange interpolation in 3D for arrays
  !============================================================================!
  function Lagrange_3d_array(x, y, z, f, point) result(solution)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(3), intent(in) :: point
    real(kind=8), dimension(:,:,:), intent(in) :: x, y, z
    real(kind=8), dimension(:,:,:,:), intent(in) :: f
    !
    ! output
    !
    real(kind=8) :: solution(size(f,1))
    !
    ! internal
    !
    integer :: i

    do i=1, size(f,1)
      solution(i) = Lagrange_3d(x,y,z,f(i,:,:,:),point)
    end do

  end function Lagrange_3d_array

end module interpolation_lib
