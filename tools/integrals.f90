module integrals
  use mpi_lib
  implicit none

  integer, private :: status(MPI_STATUS_SIZE)

  interface trapezium
    module procedure :: trapezium_real, trapezium_complex
  end interface trapezium

contains
  !============================================================================!
  real(kind=8) function trapezium_real(f, dx, dy, dz)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: dx, dy, dz
    real(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! internal
    !
    integer :: i
    real(kind=8) :: f0d, f0d_1, f0d_2
    real(kind=8), dimension(size(f,1)) :: f1d, f1d_1, f1d_2
    real(kind=8), dimension(size(f,1), size(f,2)) :: f2d

    f2d = 0.0d0
    do i=1, size(f,3)-1
      f2d = f2d + f(:,:,i) + f(:,:,i+1)
    end do

    f1d   = 0.0d0
    f1d_1 = 0.0d0
    f1d_2 = 0.0d0
    do i=1, size(f,2)-1
      f1d = f1d + f2d(:,i) + f2d(:,i+1)
    end do

    if(col_nproc>1) then
      if(col_id.eq.0) then

        call MPI_SENDRECV(&
        f2d(:,size(f,2)), size(f,1), MPI_DOUBLE_PRECISION, col_id+1, 0, &
        f1d_1           , size(f,1), MPI_DOUBLE_PRECISION, col_id+1, 1, &
        col_comm, status, ierr)

      else if(col_id.eq.col_nproc-1) then

        call MPI_SENDRECV(&
        f2d(:,1), size(f,1), MPI_DOUBLE_PRECISION, col_id-1, 1, &
        f1d_2   , size(f,1), MPI_DOUBLE_PRECISION, col_id-1, 0, &
        col_comm, status, ierr)

      else

        call MPI_SENDRECV(&
        f2d(:,size(f,2)), size(f,1), MPI_DOUBLE_PRECISION, col_id+1, 0, &
        f1d_1           , size(f,1), MPI_DOUBLE_PRECISION, col_id+1, 1, &
        col_comm, status, ierr)

        call MPI_SENDRECV(&
        f2d(:,1), size(f,1), MPI_DOUBLE_PRECISION, col_id-1, 1, &
        f1d_2   , size(f,1), MPI_DOUBLE_PRECISION, col_id-1, 0, &
        col_comm, status, ierr)

      end if
    end if

    f1d = f1d + f1d_1 + f1d_2

    f0d   = 0.0d0
    f0d_1 = 0.0d0
    f0d_2 = 0.0d0
    do i=1, size(f,1)-1
      f0d = f0d + f1d(i) + f1d(i+1)
    end do

    if(row_nproc>1) then
      if(row_id.eq.0) then

        call MPI_SENDRECV(&
        f1d(size(f,1)), 1, MPI_DOUBLE_PRECISION, row_id+1, 0, &
        f0d_1         , 1, MPI_DOUBLE_PRECISION, row_id+1, 1, &
        row_comm, status, ierr)

      else if(row_id.eq.row_nproc-1) then

        call MPI_SENDRECV(&
        f1d(1), 1, MPI_DOUBLE_PRECISION, row_id-1, 1, &
        f0d_2 , 1, MPI_DOUBLE_PRECISION, row_id-1, 0, &
        row_comm, status, ierr)

      else

        call MPI_SENDRECV(&
        f1d(size(f,1)), 1, MPI_DOUBLE_PRECISION, row_id+1, 0, &
        f0d_1         , 1, MPI_DOUBLE_PRECISION, row_id+1, 1, &
        row_comm, status, ierr)

        call MPI_SENDRECV(&
        f1d(1), 1, MPI_DOUBLE_PRECISION, row_id-1, 1, &
        f0d_2 , 1, MPI_DOUBLE_PRECISION, row_id-1, 0, &
        row_comm, status, ierr)

      end if
    end if

    f0d = 0.1250d0 * (f0d + f0d_1 + f0d_2) * dx * dy * dz

    trapezium_real = 0.0d0
    call MPI_ALLREDUCE(f0d, trapezium_real, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  end function trapezium_real
  !============================================================================!
  complex(kind=8) function trapezium_complex(f, dx, dy, dz)
    implicit none
    !
    ! input
    !
    real(kind=8), intent(in) :: dx, dy, dz
    complex(kind=8), dimension(:,:,:), intent(in) :: f
    !
    ! internal
    !
    integer :: i
    complex(kind=8) :: f0d, f0d_1, f0d_2
    complex(kind=8), dimension(size(f,1)) :: f1d, f1d_1, f1d_2
    complex(kind=8), dimension(size(f,1), size(f,2)) :: f2d

    f2d = 0.0d0
    do i=1, size(f,3)-1
      f2d = f2d + f(:,:,i) + f(:,:,i+1)
    end do

    f1d   = 0.0d0
    f1d_1 = 0.0d0
    f1d_2 = 0.0d0
    do i=1, size(f,2)-1
      f1d = f1d + f2d(:,i) + f2d(:,i+1)
    end do

    if(col_nproc>1) then
      if(col_id.eq.0) then

        call MPI_SENDRECV(&
        f2d(:,size(f,2)), size(f,1), MPI_DOUBLE_COMPLEX, col_id+1, 0, &
        f1d_1           , size(f,1), MPI_DOUBLE_COMPLEX, col_id+1, 1, &
        col_comm, status, ierr)

      else if(col_id.eq.col_nproc-1) then

        call MPI_SENDRECV(&
        f2d(:,1), size(f,1), MPI_DOUBLE_COMPLEX, col_id-1, 1, &
        f1d_2   , size(f,1), MPI_DOUBLE_COMPLEX, col_id-1, 0, &
        col_comm, status, ierr)

      else

        call MPI_SENDRECV(&
        f2d(:,size(f,2)), size(f,1), MPI_DOUBLE_COMPLEX, col_id+1, 0, &
        f1d_1           , size(f,1), MPI_DOUBLE_COMPLEX, col_id+1, 1, &
        col_comm, status, ierr)

        call MPI_SENDRECV(&
        f2d(:,1), size(f,1), MPI_DOUBLE_COMPLEX, col_id-1, 1, &
        f1d_2   , size(f,1), MPI_DOUBLE_COMPLEX, col_id-1, 0, &
        col_comm, status, ierr)

      end if
    end if

    f1d = f1d + f1d_1 + f1d_2

    f0d   = 0.0d0
    f0d_1 = 0.0d0
    f0d_2 = 0.0d0
    do i=1, size(f,1)-1
      f0d = f0d + f1d(i) + f1d(i+1)
    end do

    if(row_nproc>1) then
      if(row_id.eq.0) then

        call MPI_SENDRECV(&
        f1d(size(f,1)), 1, MPI_DOUBLE_COMPLEX, row_id+1, 0, &
        f0d_1         , 1, MPI_DOUBLE_COMPLEX, row_id+1, 1, &
        row_comm, status, ierr)

      else if(row_id.eq.row_nproc-1) then

        call MPI_SENDRECV(&
        f1d(1), 1, MPI_DOUBLE_COMPLEX, row_id-1, 1, &
        f0d_2 , 1, MPI_DOUBLE_COMPLEX, row_id-1, 0, &
        row_comm, status, ierr)

      else

        call MPI_SENDRECV(&
        f1d(size(f,1)), 1, MPI_DOUBLE_COMPLEX, row_id+1, 0, &
        f0d_1         , 1, MPI_DOUBLE_COMPLEX, row_id+1, 1, &
        row_comm, status, ierr)

        call MPI_SENDRECV(&
        f1d(1), 1, MPI_DOUBLE_COMPLEX, row_id-1, 1, &
        f0d_2 , 1, MPI_DOUBLE_COMPLEX, row_id-1, 0, &
        row_comm, status, ierr)

      end if
    end if

    f0d = 0.1250d0 * (f0d + f0d_1 + f0d_2) * dx * dy * dz

    trapezium_complex = 0.0d0
    call MPI_ALLREDUCE(f0d, trapezium_complex, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)

  end function trapezium_complex
  !============================================================================!

end module integrals
