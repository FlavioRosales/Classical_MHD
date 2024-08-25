module differential_operators
  use finite_differences
  implicit none

  interface laplacian_2
    module procedure :: laplacian_2_real, laplacian_2_complex
  end interface laplacian_2

  interface laplacian_4
    module procedure :: laplacian_4_real, laplacian_4_complex
  end interface laplacian_4

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

end module differential_operators
