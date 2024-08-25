program main
  use global_numbers
  implicit none
  !
  ! Start MPI interface
  !
  call MPI_INIT(ierr)
  !
  ! Check the number of processors used
  !
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
  !
  ! Check the global identifier of each processor
  !
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  !
  ! Create a cartesian topology and new communicators
  !
  call MPI_CREATE_CARTEIAN_COMM

  !============================================================================!
  !                     The code starts here                                   !
  !============================================================================!

  call input_data
  call initialize_variables
  call system_evolve

  call MPI_BARRIER(MPI_COMM_WORLD, ierr)


  !============================================================================!
  !                     Here ends the code                                     !
  !============================================================================!
  !
  ! End MPI interface
  !
  call MPI_FINALIZE(ierr)


end program main
