module finite_differences
  use mpi_lib, only: transpose_z_to_y, transpose_y_to_x
  use finite_differences_1d
  implicit none

  interface first_derivative_x_2
    module procedure :: first_derivative_x_2_real, first_derivative_x_2_complex
  end interface first_derivative_x_2

  interface first_derivative_x_4
    module procedure :: first_derivative_x_4_real, first_derivative_x_4_complex
  end interface first_derivative_x_4

  interface first_derivative_y_2
    module procedure :: first_derivative_y_2_real, first_derivative_y_2_complex
  end interface first_derivative_y_2

  interface first_derivative_y_4
    module procedure :: first_derivative_y_4_real, first_derivative_y_4_complex
  end interface first_derivative_y_4

  interface first_derivative_z_2
    module procedure :: first_derivative_z_2_real, first_derivative_z_2_complex
  end interface first_derivative_z_2

  interface first_derivative_z_4
    module procedure :: first_derivative_z_4_real, first_derivative_z_4_complex
  end interface first_derivative_z_4

  interface second_derivative_xx_2
    module procedure :: second_derivative_xx_2_real, second_derivative_xx_2_complex
  end interface second_derivative_xx_2

  interface second_derivative_xx_4
    module procedure :: second_derivative_xx_4_real, second_derivative_xx_4_complex
  end interface second_derivative_xx_4

  interface second_derivative_yy_2
    module procedure :: second_derivative_yy_2_real, second_derivative_yy_2_complex
  end interface second_derivative_yy_2

  interface second_derivative_yy_4
    module procedure :: second_derivative_yy_4_real, second_derivative_yy_4_complex
  end interface second_derivative_yy_4

  interface second_derivative_zz_2
    module procedure :: second_derivative_zz_2_real, second_derivative_zz_2_complex
  end interface second_derivative_zz_2

  interface second_derivative_zz_4
    module procedure :: second_derivative_zz_4_real, second_derivative_zz_4_complex
  end interface second_derivative_zz_4

  interface second_derivative_xy_2
    module procedure :: second_derivative_xy_2_real, second_derivative_xy_2_complex
  end interface second_derivative_xy_2

  interface second_derivative_xy_4
    module procedure :: second_derivative_xy_4_real, second_derivative_xy_4_complex
  end interface second_derivative_xy_4

  interface second_derivative_xz_2
    module procedure :: second_derivative_xz_2_real, second_derivative_xz_2_complex
  end interface second_derivative_xz_2

  interface second_derivative_xz_4
    module procedure :: second_derivative_xz_4_real, second_derivative_xz_4_complex
  end interface second_derivative_xz_4

contains

  !============================================================================!
  function first_derivative_x_2_real(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_real(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function first_derivative_x_2_real
  !============================================================================!
  function first_derivative_x_2_complex(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_complex(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function first_derivative_x_2_complex
  !============================================================================!
  function first_derivative_x_4_real(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_real(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function first_derivative_x_4_real
  !============================================================================!
  function first_derivative_x_4_complex(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_complex(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function first_derivative_x_4_complex
  !============================================================================!
  function first_derivative_y_2_real(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_real(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function first_derivative_y_2_real
  !============================================================================!
  function first_derivative_y_2_complex(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_complex(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function first_derivative_y_2_complex
  !============================================================================!
  function first_derivative_y_4_real(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_real(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function first_derivative_y_4_real
  !============================================================================!
  function first_derivative_y_4_complex(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_complex(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function first_derivative_y_4_complex
  !============================================================================!
  function first_derivative_z_2_real(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_real(f(i,j,:), dz)
      end do
    end do

  end function first_derivative_z_2_real
  !============================================================================!
  function first_derivative_z_2_complex(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_2_complex(f(i,j,:), dz)
      end do
    end do

  end function first_derivative_z_2_complex
  !============================================================================!
  function first_derivative_z_4_real(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_real(f(i,j,:), dz)
      end do
    end do

  end function first_derivative_z_4_real
  !============================================================================!
  function first_derivative_z_4_complex(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = first_derivative_4_complex(f(i,j,:), dz)
      end do
    end do

  end function first_derivative_z_4_complex
  !============================================================================!
  function second_derivative_xx_2_real(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_real(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function second_derivative_xx_2_real
  !============================================================================!
  function second_derivative_xx_2_complex(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_complex(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function second_derivative_xx_2_complex
  !============================================================================!
  function second_derivative_xx_4_real(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_real(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function second_derivative_xx_4_real
  !============================================================================!
  function second_derivative_xx_4_complex(f, dx) result(df)
    implicit none
    real(kind=8), intent(in) :: dx
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_y_to_x(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_4_complex(f_aux(i,j,:), dx)
      end do
    end do
    df = transpose_y_to_x(df)

  end function second_derivative_xx_4_complex
  !============================================================================!
  function second_derivative_yy_2_real(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_real(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function second_derivative_yy_2_real
  !============================================================================!
  function second_derivative_yy_2_complex(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_complex(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function second_derivative_yy_2_complex
  !============================================================================!
  function second_derivative_yy_4_real(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_4_real(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function second_derivative_yy_4_real
  !============================================================================!
  function second_derivative_yy_4_complex(f, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: f_aux, df

    integer :: i, j

    f_aux = transpose_z_to_y(f)
    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_4_complex(f_aux(i,j,:), dy)
      end do
    end do
    df = transpose_z_to_y(df)

  end function second_derivative_yy_4_complex
  !============================================================================!
  function second_derivative_zz_2_real(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_real(f(i,j,:), dz)
      end do
    end do

  end function second_derivative_zz_2_real
  !============================================================================!
  function second_derivative_zz_2_complex(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_2_complex(f(i,j,:), dz)
      end do
    end do

  end function second_derivative_zz_2_complex
  !============================================================================!
  function second_derivative_zz_4_real(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_4_real(f(i,j,:), dz)
      end do
    end do

  end function second_derivative_zz_4_real
  !============================================================================!
  function second_derivative_zz_4_complex(f, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    integer :: i, j

    do i=1, size(f,1)
      do j=1, size(f,2)
        df(i,j,:) = second_derivative_4_complex(f(i,j,:), dz)
      end do
    end do

  end function second_derivative_zz_4_complex
  !============================================================================!
  function second_derivative_xy_2_real(f, dx, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_2(f = first_derivative_y_2(f = f, dy = dy), dx = dx)

  end function second_derivative_xy_2_real
  !============================================================================!
  function second_derivative_xy_2_complex(f, dx, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_2(f = first_derivative_y_2(f = f, dy = dy), dx = dx)

  end function second_derivative_xy_2_complex
  !============================================================================!
  function second_derivative_xy_4_real(f, dx, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_4(f = first_derivative_y_4(f = f, dy = dy), dx = dx)

  end function second_derivative_xy_4_real
  !============================================================================!
  function second_derivative_xy_4_complex(f, dx, dy) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dy
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_4(f = first_derivative_y_4(f = f, dy = dy), dx = dx)

  end function second_derivative_xy_4_complex
  !============================================================================!
  function second_derivative_xz_2_real(f, dx, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_2(f = first_derivative_z_2(f = f, dz = dz), dx = dx)

  end function second_derivative_xz_2_real
  !============================================================================!
  function second_derivative_xz_2_complex(f, dx, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_2(f = first_derivative_z_2(f = f, dz = dz), dx = dx)

  end function second_derivative_xz_2_complex
  !============================================================================!
  function second_derivative_xz_4_real(f, dx, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    real(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: df

    df = first_derivative_x_4(f = first_derivative_z_4(f = f, dz = dz), dx = dx)

  end function second_derivative_xz_4_real
  !============================================================================!
  function second_derivative_xz_4_complex(f, dx, dz) result(df)
    implicit none
    real(kind=8), intent(in) :: dx, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    complex(kind=8), dimension(size(f,3), size(f,2), size(f,1)) :: df

    df = first_derivative_x_4(f = first_derivative_z_4(f = f, dz = dz), dx = dx)

  end function second_derivative_xz_4_complex
  !============================================================================!


end module finite_differences
