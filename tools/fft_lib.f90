module fft_lib
  use fft_1d_lib
  use mpi_lib
  implicit none

contains

  function sinfft(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! output
    !
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: fft_f
    !
    ! internal
    !
    integer :: i, j, k
    integer :: Nx, Ny, Nz

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3)

    !
    ! fft along the z-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(f(i,j,:))
      end do
    end do
    !
    ! tranpose z to y
    !
    fft_f = transpose_z_to_y(fft_f)
    !
    ! fft along the y-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(fft_f(i,j,:))
      end do
    end do
    !
    ! tranpose y to x
    !
    fft_f = transpose_y_to_x(fft_f)
    !
    ! fft along the x-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(fft_f(i,j,:))
      end do
    end do


  end function sinfft

  function isinfft(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! output
    !
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: fft_f
    !
    ! internal
    !
    integer :: i, j, k
    integer :: Nx, Ny, Nz

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3)

    !
    ! ifft along the x-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(f(i,j,:))
      end do
    end do
    !
    ! tranpose x to y
    !
    fft_f = transpose_y_to_x(fft_f)
    !
    ! fft along the y-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(fft_f(i,j,:))
      end do
    end do
    !
    ! tranpose y to z
    !
    fft_f = transpose_z_to_y(fft_f)
    !
    ! fft along the z-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = sinfft_1d(fft_f(i,j,:))
      end do
    end do

    fft_f = fft_f * 8.0d0 / dble(Nz)**3

  end function isinfft

  function fft(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! output
    !
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: fft_f
    !
    ! internal
    !
    integer :: i, j, k
    integer :: Nx, Ny, Nz

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3)

    !
    ! fft along the z-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fft_1d(f(i,j,:))
      end do
    end do
    !
    ! tranpose z to y
    !
    fft_f = transpose_z_to_y(fft_f)
    !
    ! fft along the y-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fft_1d(fft_f(i,j,:))
      end do
    end do
    !
    ! tranpose y to x
    !
    fft_f = transpose_y_to_x(fft_f)
    !
    ! fft along the x-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fft_1d(fft_f(i,j,:))
      end do
    end do

  end function fft

  function fftshift(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! output
    !
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: fft_f
    !
    ! internal
    !
    integer :: i, j, k
    integer :: Nx, Ny, Nz

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3)

    !
    ! fft along the z-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fftshift_1d(f(i,j,:))
      end do
    end do
    !
    ! tranpose z to y
    !
    fft_f = transpose_z_to_y(fft_f)
    !
    ! fft along the y-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fftshift_1d(fft_f(i,j,:))
      end do
    end do
    !
    ! tranpose y to x
    !
    fft_f = transpose_y_to_x(fft_f)
    !
    ! fft along the x-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = fftshift_1d(fft_f(i,j,:))
      end do
    end do

  end function fftshift

  function ifft(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! output
    !
    complex(kind=8), dimension(size(f,1), size(f,2), size(f,3)) :: fft_f
    !
    ! internal
    !
    integer :: i, j, k
    integer :: Nx, Ny, Nz

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3)

    !
    ! ifft along the x-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = ifft_1d(f(i,j,:))
      end do
    end do
    !
    ! tranpose x to y
    !
    fft_f = transpose_y_to_x(fft_f)
    !
    ! fft along the y-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = ifft_1d(fft_f(i,j,:))
      end do
    end do
    !
    ! tranpose y to z
    !
    fft_f = transpose_z_to_y(fft_f)
    !
    ! fft along the z-axis.
    !
    do i=1, Nx
      do j=1, Ny
        fft_f(i,j,:) = ifft_1d(fft_f(i,j,:))
      end do
    end do

  end function ifft

end module fft_lib
