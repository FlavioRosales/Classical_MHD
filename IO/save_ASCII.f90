module save_ASCII
  use mpi_lib
  implicit none

  logical, private :: exist
  integer, private :: i, j, k, ios
  integer, parameter, private :: iounit = 12345

  interface save_t
    module procedure :: save_t_const, save_t_vect
  end interface save_t

contains

  !============================================================================!
  subroutine save0d_const(t, f, name)
    implicit none
    real(kind=8), intent(in) :: t, f
    character(len=*), intent(in) :: name

    inquire(file = name, exist = exist)

    if(exist) then
      open(unit=iounit, file=name, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "
    else
      open(unit=iounit, file=name, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 0d"
    end if

      write(unit=iounit, fmt=*, iostat=ios) t, f
      if ( ios /= 0 ) stop "Write error in file unit "

    close(unit=iounit, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file 0d"

  end subroutine save0d_const
  !============================================================================!
  subroutine save0d_vect(t, f, name)
    implicit none
    real(kind=8), intent(in) :: t, f(:)
    character(len=*), intent(in) :: name

    inquire(file = name, exist = exist)

    if(exist) then
      open(unit=iounit, file=name, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "
    else
      open(unit=iounit, file=name, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 0d"
    end if

      write(unit=iounit, fmt=*, iostat=ios) t, f
      if ( ios /= 0 ) stop "Write error in file unit "

    close(unit=iounit, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file 0d"

  end subroutine save0d_vect
  !============================================================================!
  subroutine save1d(axis, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: axis, f
    character(len=*), intent(in) :: name

    inquire(file = name, exist = exist)

    if(exist) then
      open(unit=iounit, file=name, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "
    else
      open(unit=iounit, file=name, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 1d"
    end if

      do i=1, size(f)
        write(unit=iounit, fmt=*) axis(i), f(i)
      end do
      write(unit=iounit, fmt=*)
      write(unit=iounit, fmt=*)

    close(unit=iounit, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file 1d"

  end subroutine save1d
  !============================================================================!
  subroutine save2d(x, y, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, y, f
    character(len=*), intent(in) :: name

    inquire(file = name, exist = exist)
    if(exist) then
      open(unit=iounit, file=name, status="old", action="write", position = 'append')
    else
      open(unit=iounit, file=name, status="new", action="write")
    end if

      do i=1, size(f,1)
        do j=1, size(f,2)
          write(unit=iounit, fmt=*) x(i,j), y(i,j), f(i,j)
        end do
        write(unit=iounit, fmt=*)
      end do
      write(unit=iounit, fmt=*)
      write(unit=iounit, fmt=*)

    close(unit=iounit)

  end subroutine save2d
  !============================================================================!
  subroutine save2d_vect(axis_1, axis_2, f_1, f_2, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: axis_1, axis_2, f_1, f_2
    character(len=*), intent(in) :: name

    inquire(file = name, exist = exist)

    if(exist) then
      open(unit=iounit, file=name, iostat=ios, status="old", action="write", position = 'append')
      if ( ios /= 0 ) stop "Error opening file "
    else
      open(unit=iounit, file=name, iostat=ios, status="new", action="write")
      if ( ios /= 0 ) stop "Error opening file 2d"
    end if

      do i=1, size(f_1,1)
        do j=1, size(f_1,2)
          write(unit=iounit, fmt=*) axis_1(i,j), axis_2(i,j), f_1(i,j), f_2(i,j)
        end do
        write(unit=iounit, fmt=*)
      end do
      write(unit=iounit, fmt=*)
      write(unit=iounit, fmt=*)

    close(unit=iounit, iostat=ios)
    if ( ios /= 0 ) stop "Error closing file 2d"

  end subroutine save2d_vect
  !============================================================================!
  subroutine save_t_const(t, f, name)
    implicit none
    real(kind=8), intent(in) :: t, f
    character(len=*), intent(in) :: name

    if(rank.eq.master) then
      call save0d_const(t, f, name//'.t.asc')
    end if

  end subroutine save_t_const
  !============================================================================!
  subroutine save_t_vect(t, f, name)
    implicit none
    real(kind=8), intent(in) :: t
    real(kind=8), dimension(:), intent(in) :: f
    character(len=*), intent(in) :: name

    if(rank.eq.master) then
      call save0d_vect(t, f, name//'.t.asc')
    end if

  end subroutine save_t_vect
  !============================================================================!
  subroutine save_x(x, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: x, f
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(f)*row_nproc) :: x_global, f_global

    x_global = rebuild_x(x)
    f_global = rebuild_x(f)

    if(row_id.eq.master) then
      call save1d(x_global, f_global, name//'.x.asc')
    end if

  end subroutine save_x
  !============================================================================!
  subroutine save_y(y, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: y, f
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(f)*col_nproc) :: y_global, f_global

    y_global = rebuild_y(y)
    f_global = rebuild_y(f)

    if(col_id.eq.master) then
      call save1d(y_global, f_global, name//'.y.asc')
    end if

  end subroutine save_y
  !============================================================================!
  subroutine save_z(z, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: z, f
    character(len=*), intent(in) :: name

    call save1d(z, f, name//'.z.asc')

  end subroutine save_z
  !============================================================================!
  subroutine save_r(z, f, name)
    implicit none
    real(kind=8), dimension(:), intent(in) :: z, f
    character(len=*), intent(in) :: name

    call save1d(z, f, name//'.r.asc')

  end subroutine save_r
  !============================================================================!
  subroutine save_xy(x, y, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, y, f
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(f,1)*row_nproc, size(f,2)*col_nproc) :: x_global, y_global, f_global

    x_global = rebuild_xy(x)
    y_global = rebuild_xy(y)
    f_global = rebuild_xy(f)

    if(rank.eq.master) then
      call save2d(x_global, y_global, f_global, name // '.xy.asc')
    end if

  end subroutine save_xy
  !============================================================================!
  subroutine save_xz(x, z, f, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, z, f
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(f,1)*row_nproc, size(f,2)) :: x_global, z_global, f_global

    x_global = rebuild_xz(x)
    z_global = rebuild_xz(z)
    f_global = rebuild_xz(f)

    if(row_id.eq.master) then
      call save2d(x_global, z_global, f_global, name // '.xz.asc')
    end if

  end subroutine save_xz
  !============================================================================!
  subroutine save_xy_vect(x, y, fx, fy, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, y, fx, fy
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(fx,1)*row_nproc, size(fx,2)*col_nproc) :: x_global, y_global, fx_global, fy_global

    x_global = rebuild_xy(x)
    y_global = rebuild_xy(y)
    fx_global = rebuild_xy(fx)
    fy_global = rebuild_xy(fy)

    if(rank.eq.master) then
      call save2d_vect(x_global, y_global, fx_global, fy_global, name // '.xy.asc')
    end if

  end subroutine save_xy_vect
  !============================================================================!
  subroutine save_xz_vect(x, z, fx, fz, name)
    implicit none
    real(kind=8), dimension(:,:), intent(in) :: x, z, fx, fz
    character(len=*), intent(in) :: name

    real(kind=8), dimension(size(fx,1)*row_nproc, size(fx,2)) :: x_global, z_global, fx_global, fz_global

    x_global = rebuild_xz(x)
    z_global = rebuild_xz(z)
    fx_global = rebuild_xz(fx)
    fz_global = rebuild_xz(fz)

    if(row_id.eq.master) then
      call save2d_vect(x_global, z_global, fx_global, fz_global, name // '.xz.asc')
    end if

  end subroutine save_xz_vect
  !============================================================================!

end module save_ASCII
