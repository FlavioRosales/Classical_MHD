module periodic_finite_differences_1d
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

    df(0) = (- f(N) + f(1)) * res
    do i=1, N-1
      df(i) = (- f(i-1) + f(i+1)) * res
    end do
    df(N) = (- f(N-1) + f(0)) * res

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

    df(0) = (- f(N) + f(1)) * res
    do i=1, N-1
      df(i) = (- f(i-1) + f(i+1)) * res
    end do
    df(N) = (- f(N-1) + f(0)) * res

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

    df(0) = (f(N-1) - 8.0d0 * f(N) + 8.0d0 * f(1) - f(2)) * res
    df(1) = (f(N-0) - 8.0d0 * f(0) + 8.0d0 * f(2) - f(3)) * res
    do i=2, N-2
      df(i) = (f(i-2) - 8.0d0 * f(i-1) + 8.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (f(N-3) - 8.0d0 * f(N-2) + 8.0d0 * f(N) - f(0)) * res
    df(N-0) = (f(N-2) - 8.0d0 * f(N-1) + 8.0d0 * f(0) - f(1)) * res

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

    df(0) = (f(N-1) - 8.0d0 * f(N) + 8.0d0 * f(1) - f(2)) * res
    df(1) = (f(N-0) - 8.0d0 * f(0) + 8.0d0 * f(2) - f(3)) * res
    do i=2, N-2
      df(i) = (f(i-2) - 8.0d0 * f(i-1) + 8.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (f(N-3) - 8.0d0 * f(N-2) + 8.0d0 * f(N) - f(0)) * res
    df(N-0) = (f(N-2) - 8.0d0 * f(N-1) + 8.0d0 * f(0) - f(1)) * res


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

    df(0) = (f(N) - 2.0d0 * f(0) + f(1)) * res
    do i=1, N-1
      df(i) = (f(i-1) - 2.0d0 * f(i) + f(i+1)) * res
    end do
    df(N) = (f(N-1) - 2.0d0 * f(N) + f(0)) * res

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

    df(0) = (f(N) - 2.0d0 * f(0) + f(1)) * res
    do i=1, N-1
      df(i) = (f(i-1) - 2.0d0 * f(i) + f(i+1)) * res
    end do
    df(N) = (f(N-1) - 2.0d0 * f(N) + f(0)) * res

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

    df(0) = (-f(N-1) + 16.0d0 * f(N) - 30.0d0 * f(0) + 16.0d0 * f(1) - f(2)) * res
    df(1) = (-f(N-0) + 16.0d0 * f(0) - 30.0d0 * f(1) + 16.0d0 * f(2) - f(3)) * res
    do i=2, N-2
      df(i) = (-f(i-2) + 16.0d0 * f(i-1) - 30.0d0 * f(i) + 16.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (-f(N-3) + 16.0d0 * f(N-2) - 30.0d0 * f(N-1) + 16.0d0 * f(N) - f(0)) * res
    df(N-0) = (-f(N-2) + 16.0d0 * f(N-1) - 30.0d0 * f(N-0) + 16.0d0 * f(0) - f(1)) * res

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

    df(0) = (-f(N-1) + 16.0d0 * f(N) - 30.0d0 * f(0) + 16.0d0 * f(1) - f(2)) * res
    df(1) = (-f(N-0) + 16.0d0 * f(0) - 30.0d0 * f(1) + 16.0d0 * f(2) - f(3)) * res
    do i=2, N-2
      df(i) = (-f(i-2) + 16.0d0 * f(i-1) - 30.0d0 * f(i) + 16.0d0 * f(i+1) - f(i+2)) * res
    end do
    df(N-1) = (-f(N-3) + 16.0d0 * f(N-2) - 30.0d0 * f(N-1) + 16.0d0 * f(N) - f(0)) * res
    df(N-0) = (-f(N-2) + 16.0d0 * f(N-1) - 30.0d0 * f(N-0) + 16.0d0 * f(0) - f(1)) * res

  end function second_derivative_4_complex
  !============================================================================!
  function second_derivative_8_real(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    real(kind=8), dimension(0:), intent(in) :: f
    real(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (5040.0d0 * h**2)

    df(0) = (-9.0d0*f(N-3)+128.0d0*f(N-2)-1008.0d0*f(N-1)+8064.0d0*f(N)-&
    14350.0d0*f(0)+8064.0d0*f(1)-1008.0d0*f(2)+128.0d0*f(3)-9.0d0*f(4)) * res
    df(1) = (-9.0d0*f(N-2)+128.0d0*f(N-1)-1008.0d0*f(N-0)+8064.0d0*f(0)-&
    14350.0d0*f(1)+8064.0d0*f(2)-1008.0d0*f(3)+128.0d0*f(4)-9.0d0*f(5)) * res
    df(2) = (-9.0d0*f(N-1)+128.0d0*f( N )-1008.0d0*f( 0 )+8064.0d0*f(1)-&
    14350.0d0*f(2)+8064.0d0*f(3)-1008.0d0*f(4)+128.0d0*f(5)-9.0d0*f(6)) * res
    df(3) = (-9.0d0*f(N-0)+128.0d0*f( 0 )-1008.0d0*f( 1 )+8064.0d0*f(2)-&
    14350.0d0*f(3)+8064.0d0*f(4)-1008.0d0*f(5)+128.0d0*f(6)-9.0d0*f(7)) * res
    do i=4, N-4
      df(i) = (-9.0d0*f(i-4)+128.0d0*f(i-3)-1008.0d0*f(i-2)+8064.0d0*f(i-1)-&
      14350.0d0*f(i)+8064.0d0*f(i+1)-1008.0d0*f(i+2)+128.0d0*f(i+3)-9.0d0*f(i+4)) * res
    end do
    df(N-3) = (-9.0d0*f(N-7)+128.0d0*f(N-6)-1008.0d0*f(N-5)+8064.0d0*f(N-4)-&
    14350.0d0*f(N-3)+8064.0d0*f(N-2)-1008.0d0*f(N-1)+128.0d0*f(N)-9.0d0*f(0)) * res
    df(N-2) = (-9.0d0*f(N-6)+128.0d0*f(N-5)-1008.0d0*f(N-4)+8064.0d0*f(N-3)-&
    14350.0d0*f(N-2)+8064.0d0*f(N-1)-1008.0d0*f( N )+128.0d0*f(0)-9.0d0*f(1)) * res
    df(N-1) = (-9.0d0*f(N-5)+128.0d0*f(N-4)-1008.0d0*f(N-3)+8064.0d0*f(N-2)-&
    14350.0d0*f(N-1)+8064.0d0*f( N )-1008.0d0*f( 0 )+128.0d0*f(1)-9.0d0*f(2)) * res
    df( N ) = (-9.0d0*f(N-4)+128.0d0*f(N-3)-1008.0d0*f(N-2)+8064.0d0*f(N-1)-&
    14350.0d0*f( N )+8064.0d0*f( 0 )-1008.0d0*f( 1 )+128.0d0*f(2)-9.0d0*f(3)) * res

  end function second_derivative_8_real
  !============================================================================!
  function second_derivative_8_complex(f, h) result(df)
    implicit none
    real(kind=8), intent(in) :: h
    complex(kind=8), dimension(0:), intent(in) :: f
    complex(kind=8), dimension(0:size(f)-1) :: df

    integer :: i, N
    real(kind=8) :: res

    N = size(f)-1; res = 1.0d0 / (5040.0d0 * h**2)

    df(0) = (-9.0d0*f(N-3)+128.0d0*f(N-2)-1008.0d0*f(N-1)+8064.0d0*f(N)-&
    14350.0d0*f(0)+8064.0d0*f(1)-1008.0d0*f(2)+128.0d0*f(3)-9.0d0*f(4)) * res
    df(1) = (-9.0d0*f(N-2)+128.0d0*f(N-1)-1008.0d0*f(N-0)+8064.0d0*f(0)-&
    14350.0d0*f(1)+8064.0d0*f(2)-1008.0d0*f(3)+128.0d0*f(4)-9.0d0*f(5)) * res
    df(2) = (-9.0d0*f(N-1)+128.0d0*f( N )-1008.0d0*f( 0 )+8064.0d0*f(1)-&
    14350.0d0*f(2)+8064.0d0*f(3)-1008.0d0*f(4)+128.0d0*f(5)-9.0d0*f(6)) * res
    df(3) = (-9.0d0*f(N-0)+128.0d0*f( 0 )-1008.0d0*f( 1 )+8064.0d0*f(2)-&
    14350.0d0*f(3)+8064.0d0*f(4)-1008.0d0*f(5)+128.0d0*f(6)-9.0d0*f(7)) * res
    do i=4, N-4
      df(i) = (-9.0d0*f(i-4)+128.0d0*f(i-3)-1008.0d0*f(i-2)+8064.0d0*f(i-1)-&
      14350.0d0*f(i)+8064.0d0*f(i+1)-1008.0d0*f(i+2)+128.0d0*f(i+3)-9.0d0*f(i+4)) * res
    end do
    df(N-3) = (-9.0d0*f(N-7)+128.0d0*f(N-6)-1008.0d0*f(N-5)+8064.0d0*f(N-4)-&
    14350.0d0*f(N-3)+8064.0d0*f(N-2)-1008.0d0*f(N-1)+128.0d0*f(N)-9.0d0*f(0)) * res
    df(N-2) = (-9.0d0*f(N-6)+128.0d0*f(N-5)-1008.0d0*f(N-4)+8064.0d0*f(N-3)-&
    14350.0d0*f(N-2)+8064.0d0*f(N-1)-1008.0d0*f( N )+128.0d0*f(0)-9.0d0*f(1)) * res
    df(N-1) = (-9.0d0*f(N-5)+128.0d0*f(N-4)-1008.0d0*f(N-3)+8064.0d0*f(N-2)-&
    14350.0d0*f(N-1)+8064.0d0*f( N )-1008.0d0*f( 0 )+128.0d0*f(1)-9.0d0*f(2)) * res
    df( N ) = (-9.0d0*f(N-4)+128.0d0*f(N-3)-1008.0d0*f(N-2)+8064.0d0*f(N-1)-&
    14350.0d0*f( N )+8064.0d0*f( 0 )-1008.0d0*f( 1 )+128.0d0*f(2)-9.0d0*f(3)) * res

  end function second_derivative_8_complex
  !============================================================================!

end module periodic_finite_differences_1d
