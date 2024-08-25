module FFT_1D_lib
  implicit none

  real(kind=8), parameter, private :: pi = acos(-1.0d0)
  complex(kind=8), parameter, private :: i_unreal = (0.0d0, 1.0d0)

contains

  function DFT_1D(f) result(DFT_f)
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: DFT_f
    !
    ! internal
    !
    integer :: i, j, N
    complex(kind=8) :: W

    N = size(f)
    W = exp(2.0d0 * pi * i_unreal / dble(N))

    DFT_f = 0.0d0
    do i=0, N-1
      do j=0, N-1
        DFT_f(i) = DFT_f(i) + W**(i * j) * f(j)
      end do
    end do

  end function DFT_1D

  recursive function FFT_1D(f) result(FFT_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: FFT_f
    !
    ! internal
    !
    integer :: i, N
    complex(kind=8) :: W
    complex(kind=8), dimension(0:size(f)/2-1) :: FFT_f_even, FFT_f_odd

    N = size(f)
    W = exp(2.0d0 * pi * i_unreal / dble(N))

    if(mod(N,2).eq.0) then

      FFT_f_even = FFT_1D(f = f(0:N-1:2))
      FFT_f_odd  = FFT_1D(f = f(1:N-1:2))

      do i=0, N/2-1
        FFT_f(i)       = FFT_f_even(i) + W**i * FFT_f_odd(i)
        FFT_f(i + N/2) = FFT_f_even(i) - W**i * FFT_f_odd(i)
      end do

    else

      FFT_f = DFT_1D(f)

    end if

  end function FFT_1D

  recursive function IFFT_1D(f) result(FFT_f)

    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: FFT_f

    FFT_f = conjg(FFT_1d(conjg(f))) / dble(size(f))

  end function IFFT_1D

  function fftshift_1d(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: fft_f
    !
    ! internal
    !
    integer :: N
    complex(kind=8), dimension(0:size(f)-1) :: fft

    N = size(f)/2
    fft = FFT_1D(f)

    fft_f(0:N-1) = fft(N:2*N-1)
    fft_f(N:2*N-1) = fft(0:N-1)

  end function fftshift_1d

  function ifftshift_1d(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: fft_f
    !
    ! internal
    !
    integer :: N
    complex(kind=8), dimension(0:size(f)-1) :: fft

    N = size(f)/2

    fft(0:N-1) = f(N:2*N-1)
    fft(N:2*N-1) = f(0:N-1)

    fft = FFT_1D(fft)

    fft_f(0:N-1) = fft(N:2*N-1)
    fft_f(N:2*N-1) = fft(0:N-1)

  end function ifftshift_1d

  function sinfft_1d(f) result(fft_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f)-1) :: fft_f
    !
    ! internal
    !
    integer :: j, k, N
    complex(kind=8), dimension(0:size(f)-1) :: f_temp, fft_temp

    N = size(f)

    f_temp(0) = 0.0d0
    do j=1, N-1
      f_temp(j) = sin(dble(j)*pi/dble(N)) * (f(j) + f(N-j)) + 0.50d0 * (f(j) - f(N-j))
    end do

    fft_temp = fft_1d(f_temp)

    fft_f(0) = dimag(fft_temp(0))
    fft_f(1) = 0.50d0 * dreal(fft_temp(0))

    do k=1, N/2-1
      fft_f(2*k) = dimag(fft_temp(k))
      fft_f(2*k+1) = fft_f(2*k-1) + dreal(fft_temp(k))
    end do

  end function sinfft_1d

end module FFT_1D_lib
