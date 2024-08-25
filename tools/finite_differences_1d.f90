module finite_differences_1d
  implicit none

contains
  !============================================================================!
  function first_derivative_2_real(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    real(kind=8), dimension(0:), intent(in) :: f
    real(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 0.50d0 / h

    df(0) = (-3.0d0 * f(0) + 4.0d0 * f(1) - f(2)) * res
    do i=1, N-1
      df(i) = (f(i+1) - f(i-1)) * res
    end do
    df(N) = (+3.0d0 * f(N) - 4.0d0 * f(N-1) + f(N-2)) * res

  end function first_derivative_2_real
  !============================================================================!
  function first_derivative_2_complex(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    complex(kind=8), dimension(0:), intent(in) :: f
    complex(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 0.50d0 / h

    df(0) = (-3.0d0 * f(0) + 4.0d0 * f(1) - f(2)) * res
    do i=1, N-1
      df(i) = (f(i+1) - f(i-1)) * res
    end do
    df(N) = (+3.0d0 * f(N) - 4.0d0 * f(N-1) + f(N-2)) * res

  end function first_derivative_2_complex
  !============================================================================!
  function first_derivative_4_real(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    real(kind=8), dimension(0:), intent(in) :: f
    real(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (12.0d0 * h)

    df(0) = (-25.0d0 * f(0) + 48.0d0 * f(1) - &
    36.0d0 * f(2) + 16.0d0 * f(3) - 3.0d0 * f(4)) * res
    df(1) = (-25.0d0 * f(1) + 48.0d0 * f(2) - &
    36.0d0 * f(3) + 16.0d0 * f(4) - 3.0d0 * f(5)) * res
    do i=2, N-2
      df(i) = (f(i-2) - 8.0d0 * f(i-1) + 8.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (+25.0d0 * f(N-1) - 48.0d0 * f(N-2) + &
    36.0d0 * f(N-3) - 16.0d0 * f(N-4) + 3.0d0 * f(N-5)) * res
    df(N-0) = (+25.0d0 * f(N-0) - 48.0d0 * f(N-1) + &
    36.0d0 * f(N-2) - 16.0d0 * f(N-3) + 3.0d0 * f(N-4)) * res


  end function first_derivative_4_real
  !============================================================================!
  function first_derivative_4_complex(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    complex(kind=8), dimension(0:), intent(in) :: f
    complex(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (12.0d0 * h)

    df(0) = (-25.0d0 * f(0) + 48.0d0 * f(1) - &
    36.0d0 * f(2) + 16.0d0 * f(3) - 3.0d0 * f(4)) * res
    df(1) = (-25.0d0 * f(1) + 48.0d0 * f(2) - &
    36.0d0 * f(3) + 16.0d0 * f(4) - 3.0d0 * f(5)) * res
    do i=2, N-2
      df(i) = (f(i-2) - 8.0d0 * f(i-1) + 8.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (+25.0d0 * f(N-1) - 48.0d0 * f(N-2) + &
    36.0d0 * f(N-3) - 16.0d0 * f(N-4) + 3.0d0 * f(N-5)) * res
    df(N-0) = (+25.0d0 * f(N-0) - 48.0d0 * f(N-1) + &
    36.0d0 * f(N-2) - 16.0d0 * f(N-3) + 3.0d0 * f(N-4)) * res


  end function first_derivative_4_complex
  !============================================================================!
  function second_derivative_2_real(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    real(kind=8), dimension(0:), intent(in) :: f
    real(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / h**2

    df(0) = (2.0d0 * f(0) - 5.0d0 * f(1) + 4.0d0 * f(2) - f(3)) * res
    do i=1, N-1
      df(i) = (f(i-1) - 2.0d0 * f(i) + f(i+1)) * res
    end do
    df(N) = (2.0d0 * f(N) - 5.0d0 * f(N-1) + 4.0d0 * f(N-2) - f(N-3)) * res

  end function second_derivative_2_real
  !============================================================================!
  function second_derivative_2_complex(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    complex(kind=8), dimension(0:), intent(in) :: f
    complex(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / h**2

    df(0) = (2.0d0 * f(0) - 5.0d0 * f(1) + 4.0d0 * f(2) - f(3)) * res
    do i=1, N-1
      df(i) = (f(i+1) - 2.0d0 * f(i) + f(i-1)) * res
    end do
    df(N) = (2.0d0 * f(N) - 5.0d0 * f(N-1) + 4.0d0 * f(N-2) - f(N-3)) * res

  end function second_derivative_2_complex
  !============================================================================!
  function second_derivative_4_real(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    real(kind=8), dimension(0:), intent(in) :: f
    real(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (12.0d0 * h**2)

    df(0) = (35.0d0*f(0)-104.0d0*f(1)+114.0d0*f(2)-56.0d0*f(3)+11.0d0*f(4)) * res
    df(1) = (11.0d0*f(0)-20.0d0*f(1)+6.0d0*f(2)+4.0d0*f(3)-f(4)) * res
    do i=2, N-2
      df(i) = (-f(i-2) + 16.0d0 * f(i-1) - 30.0d0 * f(i) + 16.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (11.0d0*f(N)-20.0d0*f(N-1)+6.0d0*f(N-2)+4.0d0*f(N-3)-f(N-4)) * res
    df(N) = (35.0d0*f(N)-104.0d0*f(N-1)+114.0d0*f(N-2)-56.0d0*f(N-3)+11.0d0*f(N-4)) * res

  end function second_derivative_4_real
  !============================================================================!
  function second_derivative_4_complex(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    complex(kind=8), dimension(0:), intent(in) :: f
    complex(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (12.0d0 * h**2)

    df(0) = (35.0d0*f(0)-104.0d0*f(1)+114.0d0*f(2)-56.0d0*f(3)+11.0d0*f(4)) * res
    df(1) = (11.0d0*f(0)-20.0d0*f(1)+6.0d0*f(2)+4.0d0*f(3)-f(4)) * res
    do i=2, N-2
      df(i) = (-f(i-2) + 16.0d0 * f(i-1) - 30.0d0 * f(i) + 16.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (11.0d0*f(N)-20.0d0*f(N-1)+6.0d0*f(N-2)+4.0d0*f(N-3)-f(N-4)) * res
    df(N) = (35.0d0*f(N)-104.0d0*f(N-1)+114.0d0*f(N-2)-56.0d0*f(N-3)+11.0d0*f(N-4)) * res

  end function second_derivative_4_complex
  !============================================================================!

end module finite_differences_1d
