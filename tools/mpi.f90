module mpi_lib
  use mpi
  implicit none

  integer :: ierr
  integer :: nproc
  integer :: rank
  integer :: COMM_CART
  integer, parameter :: ndims = 2
  integer, dimension(2) :: dims
  integer :: row_nproc = 1, col_nproc = 1
  logical, dimension(2) :: notperiodic = [.false., .false.]
  integer :: master = 0
  integer, dimension(2) :: coords

  !
  ! new communicators
  !
  integer :: row_comm, col_comm
  integer :: row_id, col_id
  integer, dimension(2) :: row_coords, col_coords

  interface transpose_z_to_y
    module procedure :: transpose_z_to_y_real, transpose_z_to_y_complex
  end interface transpose_z_to_y

  interface transpose_y_to_x
    module procedure :: transpose_y_to_x_real, transpose_y_to_x_complex
  end interface transpose_y_to_x

contains

  !============================================================================!
  subroutine MPI_CREATE_CARTEIAN_COMM
    implicit none
    integer :: i
    integer, allocatable, dimension(:) :: remdims
    !
    ! Creates a new communicator with cartesian topology.
    !
    call MPI_DIMS_CREATE(nproc, ndims, dims, ierr)
    call MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, notperiodic, .false., COMM_CART, ierr)
    !
    ! converts process rank to proces grid coordinates.
    !
    call MPI_CART_COORDS(COMM_CART, rank, ndims, coords, ierr)

    row_nproc = dims(1)
    col_nproc = dims(2)

    !
    ! Create a row communicator
    !
    call MPI_CART_SUB(COMM_CART, [.true., .false.], row_comm, ierr)
    call MPI_COMM_RANK(row_comm, row_id, ierr)
    call MPI_CART_COORDS(row_comm, row_id, 1, row_coords, ierr)
    !
    ! Create a column communicator
    !
    call MPI_CART_SUB(COMM_CART, [.false., .true.], col_comm, ierr)
    call MPI_COMM_RANK(col_comm, col_id, ierr)
    call MPI_CART_COORDS(col_comm, col_id, 1, col_coords, ierr)

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  end subroutine MPI_CREATE_CARTEIAN_COMM
  !============================================================================!
  function convert_rank_to_coords(local_rank) result(local_coords)
    implicit none
    integer, intent(in) :: local_rank
    integer, dimension(2) :: local_coords

    call MPI_CART_COORDS(COMM_CART, local_rank, 2, local_coords, ierr)

  end function convert_rank_to_coords
  !============================================================================!
  function convert_coords_to_rank(local_coords) result(local_rank)
    implicit none
    integer, dimension(2), intent(in) :: local_coords
    integer :: local_rank

    call MPI_CART_RANK(COMM_CART, local_coords, local_rank, ierr)

  end function convert_coords_to_rank
  !============================================================================!
  function transpose_z_to_y_complex(f) result(transpose_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:,0:,0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: transpose_f
    !
    ! internal
    !
    integer :: i, j, k, iblock, blocks
    integer :: Nx, Ny, Nz
    complex(kind=8), dimension(0:size(f,2)*size(f,3)-1) :: Ctmp, Dtmp

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3); blocks = col_nproc

    do i=0, Nx-1
      !
      ! Divide f into blocks
      !
      do j=0, Ny-1
        do iblock=0, blocks-1
          do k=0, Ny-1
            Ctmp(iblock*Ny*Ny+j*Ny+k) = f(i,j,iblock*Ny+k)
          end do
        end do
      end do
      !
      ! local transposition
      !
      do iblock=0, blocks-1
        do j=0, Ny-1
          do k=0, Ny-1
            Dtmp(iblock*Ny*Ny+j*Ny+k) = Ctmp(iblock*Ny*Ny+k*Ny+j)
          end do
        end do
      end do
      !
      ! All to all comm
      !
      call MPI_ALLTOALL(&
      Dtmp, Ny*Ny, MPI_DOUBLE_COMPLEX, &
      Ctmp, Ny*Ny, MPI_DOUBLE_COMPLEX, &
      col_comm, ierr)
      !
      ! merge blocks
      !
      do j=0, Ny-1
        do iblock=0, blocks-1
          do k=0, Ny-1
            transpose_f(i,j,iblock*Ny+k) = Ctmp(iblock*Ny*Ny+j*Ny+k)
          end do
        end do
      end do
    end do

  end function transpose_z_to_y_complex
  !============================================================================!
  function transpose_y_to_x_complex(f) result(transpose_f)
    implicit none
    !
    ! input
    !
    complex(kind=8), dimension(0:,0:,0:) :: f
    !
    ! output
    !
    complex(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: transpose_f
    !
    ! internal
    !
    integer :: i, j, k, iblock, blocks
    integer :: Nx, Ny, Nz
    complex(kind=8), dimension(0:size(f,2)*size(f,3)-1) :: Ctmp, Dtmp

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3); blocks = row_nproc

    do i=0, Ny-1
      !
      ! Divide f into blocks
      !
      do j=0, Nx-1
        do iblock=0, blocks-1
          do k=0, Nx-1
            Ctmp(iblock*Nx*Nx+j*Nx+k) = f(j,i,iblock*Nx+k)
          end do
        end do
      end do
      !
      ! local transposition
      !
      do iblock=0, blocks-1
        do j=0, Nx-1
          do k=0, Nx-1
            Dtmp(iblock*Nx*Nx+j*Nx+k) = Ctmp(iblock*Nx*Nx+k*Nx+j)
          end do
        end do
      end do
      !
      ! All to all comm
      !
      call MPI_ALLTOALL(&
      Dtmp, Nx*Nx, MPI_DOUBLE_COMPLEX, &
      Ctmp, Nx*Nx, MPI_DOUBLE_COMPLEX, &
      row_comm, ierr)
      !
      ! merge blocks
      !
      do j=0, Nx-1
        do iblock=0, blocks-1
          do k=0, Nx-1
            transpose_f(j,i,iblock*Nx+k) = Ctmp(iblock*Nx*Nx+j*Nx+k)
          end do
        end do
      end do
    end do

  end function transpose_y_to_x_complex
  !============================================================================!
  function transpose_z_to_y_real(f) result(transpose_f)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(0:,0:,0:) :: f
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: transpose_f
    !
    ! internal
    !
    integer :: i, j, k, iblock, blocks
    integer :: Nx, Ny, Nz
    real(kind=8), dimension(0:size(f,2)*size(f,3)-1) :: Ctmp, Dtmp

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3); blocks = col_nproc

    do i=0, Nx-1
      !
      ! Divide f into blocks
      !
      do j=0, Ny-1
        do iblock=0, blocks-1
          do k=0, Ny-1
            Ctmp(iblock*Ny*Ny+j*Ny+k) = f(i,j,iblock*Ny+k)
          end do
        end do
      end do
      !
      ! local transposition
      !
      do iblock=0, blocks-1
        do j=0, Ny-1
          do k=0, Ny-1
            Dtmp(iblock*Ny*Ny+j*Ny+k) = Ctmp(iblock*Ny*Ny+k*Ny+j)
          end do
        end do
      end do
      !
      ! All to all comm
      !
      call MPI_ALLTOALL(&
      Dtmp, Ny*Ny, MPI_DOUBLE_PRECISION, &
      Ctmp, Ny*Ny, MPI_DOUBLE_PRECISION, &
      col_comm, ierr)
      !
      ! merge blocks
      !
      do j=0, Ny-1
        do iblock=0, blocks-1
          do k=0, Ny-1
            transpose_f(i,j,iblock*Ny+k) = Ctmp(iblock*Ny*Ny+j*Ny+k)
          end do
        end do
      end do
    end do

  end function transpose_z_to_y_real
  !============================================================================!
  function transpose_y_to_x_real(f) result(transpose_f)
    implicit none
    !
    ! input
    !
    real(kind=8), dimension(0:,0:,0:) :: f
    !
    ! output
    !
    real(kind=8), dimension(0:size(f,1)-1, 0:size(f,2)-1, 0:size(f,3)-1) :: transpose_f
    !
    ! internal
    !
    integer :: i, j, k, iblock, blocks
    integer :: Nx, Ny, Nz
    real(kind=8), dimension(0:size(f,2)*size(f,3)-1) :: Ctmp, Dtmp

    Nx = size(f,1); Ny = size(f,2); Nz = size(f,3); blocks = row_nproc

    do i=0, Ny-1
      !
      ! Divide f into blocks
      !
      do j=0, Nx-1
        do iblock=0, blocks-1
          do k=0, Nx-1
            Ctmp(iblock*Nx*Nx+j*Nx+k) = f(j,i,iblock*Nx+k)
          end do
        end do
      end do
      !
      ! local transposition
      !
      do iblock=0, blocks-1
        do j=0, Nx-1
          do k=0, Nx-1
            Dtmp(iblock*Nx*Nx+j*Nx+k) = Ctmp(iblock*Nx*Nx+k*Nx+j)
          end do
        end do
      end do
      !
      ! All to all comm
      !
      call MPI_ALLTOALL(&
      Dtmp, Nx*Nx, MPI_DOUBLE_PRECISION, &
      Ctmp, Nx*Nx, MPI_DOUBLE_PRECISION, &
      row_comm, ierr)
      !
      ! merge blocks
      !
      do j=0, Nx-1
        do iblock=0, blocks-1
          do k=0, Nx-1
            transpose_f(j,i,iblock*Nx+k) = Ctmp(iblock*Nx*Nx+j*Nx+k)
          end do
        end do
      end do
    end do

  end function transpose_y_to_x_real
  !============================================================================!
  function rebuild_x(f) result(f_global)
    implicit none
    real(kind=8), dimension(:) :: f
    real(kind=8), dimension(size(f)*row_nproc) :: f_global
    real(kind=8), dimension(row_nproc, size(f)) :: f_local
    integer :: i
    integer :: status(MPI_STATUS_SIZE)

    if(row_id.ne.master) call MPI_SEND(f, size(f), MPI_DOUBLE_PRECISION, master, row_id, row_comm, ierr)

    call MPI_BARRIER(row_comm, ierr)

    if(row_id.eq.master) then
      f_global(1:size(f)) = f
      do i=1, row_nproc-1
        call MPI_RECV(f_global(size(f,1) * i + 1:size(f,1) * (i+1)), &
        size(f), MPI_DOUBLE_PRECISION, i, i, row_comm, status, ierr)
      end do
    end if

    call MPI_BARRIER(row_comm, ierr)

  end function rebuild_x
  !============================================================================!
  function rebuild_y(f) result(f_global)
    implicit none
    real(kind=8), dimension(:) :: f
    real(kind=8), dimension(size(f)*col_nproc) :: f_global
    real(kind=8), dimension(col_nproc, size(f)) :: f_local
    integer :: i
    integer :: status(MPI_STATUS_SIZE)

    if(col_id.ne.master) call MPI_SEND(f, size(f), MPI_DOUBLE_PRECISION, master, col_id, col_comm, ierr)

    call MPI_BARRIER(col_comm, ierr)

    if(col_id.eq.master) then
      f_global(1:size(f)) = f
      do i=1, col_nproc-1
        call MPI_RECV(f_global(size(f,1) * i + 1:size(f,1) * (i+1)), &
        size(f), MPI_DOUBLE_PRECISION, i, i, col_comm, status, ierr)
      end do
    end if

    call MPI_BARRIER(col_comm, ierr)

  end function rebuild_y
  !============================================================================!
  function rebuild_xy(f) result(f_global)
    implicit none
    real(kind=8), dimension(:,:) :: f
    real(kind=8), dimension(size(f,1)*row_nproc, size(f,2)*col_nproc) :: f_global
    real(kind=8), dimension(nproc, size(f,1), size(f,2)) :: f_local
    integer :: i, coords_i(2)
    integer :: status(MPI_STATUS_SIZE)

    f_global = 0.0d0

    if(rank.ne.master) call MPI_SEND(f, size(f), MPI_DOUBLE_PRECISION, master, rank, MPI_COMM_WORLD, ierr)

    if(rank.eq.master) then
      f_global(1:size(f,1), 1:size(f,2)) = f
      do i=1, nproc-1
        coords_i = convert_rank_to_coords(i)
        call MPI_RECV(f_global(&
        size(f,1) * coords_i(1) + 1:size(f,1) * (coords_i(1)+1), &
        size(f,2) * coords_i(2) + 1:size(f,2) * (coords_i(2)+1)), size(f), MPI_DOUBLE_PRECISION, i, i, MPI_COMM_WORLD, status, ierr)
      end do
    end if

  end function rebuild_xy
  !============================================================================!
  function rebuild_xz(f) result(f_global)
    implicit none
    real(kind=8), dimension(:,:) :: f
    real(kind=8), dimension(size(f,1)*row_nproc, size(f,2)) :: f_global
    real(kind=8), dimension(row_nproc, size(f,1), size(f,2)) :: f_local
    integer :: i
    integer :: status(MPI_STATUS_SIZE)

    if(row_id.ne.master) call MPI_SEND(f, size(f), MPI_DOUBLE_PRECISION, master, row_id, row_comm, ierr)

    if(row_id.eq.master) then
      f_global(1:size(f,1), 1:size(f,2)) = f
      do i=1, row_nproc-1
        call MPI_RECV(f_global(size(f,1) * i + 1:size(f,1) * (i+1), :), &
        size(f), MPI_DOUBLE_PRECISION, i, i, row_comm, status, ierr)
      end do
    end if

  end function rebuild_xz
  !============================================================================!

end module mpi_lib
