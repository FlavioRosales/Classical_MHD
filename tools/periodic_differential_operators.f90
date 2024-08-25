module periodic_differential_operators
  use periodic_finite_differences
  implicit none

  interface laplacian_2
    module procedure :: laplacian_2_real, laplacian_2_complex
  end interface laplacian_2

  interface laplacian_4
    module procedure :: laplacian_4_real, laplacian_4_complex
  end interface laplacian_4

  interface laplacian_8
    module procedure :: laplacian_8_real, laplacian_8_complex
  end interface laplacian_8

contains
  !============================================================================!
  function laplacian_2_real(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_2(f, dx) + second_derivative_yy_2(f, dy) + second_derivative_zz_2(f, dz)

  end function laplacian_2_real
  !============================================================================!
  function laplacian_2_complex(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_2(f, dx) + second_derivative_yy_2(f, dy) + second_derivative_zz_2(f, dz)

  end function laplacian_2_complex
  !============================================================================!
  function laplacian_4_real(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_4(f, dx) + second_derivative_yy_4(f, dy) + second_derivative_zz_4(f, dz)

  end function laplacian_4_real
  !============================================================================!
  function laplacian_4_complex(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_4(f, dx) + second_derivative_yy_4(f, dy) + second_derivative_zz_4(f, dz)

  end function laplacian_4_complex
  !============================================================================!
  function laplacian_8_real(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_4(f, dx) + second_derivative_yy_8(f, dy) + second_derivative_zz_8(f, dz)

  end function laplacian_8_real
  !============================================================================!
  function laplacian_8_complex(f, dx, dy, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = second_derivative_xx_8(f, dx) + second_derivative_yy_8(f, dy) + second_derivative_zz_8(f, dz)

  end function laplacian_8_complex
  !============================================================================!

end module periodic_differential_operators
