module PopulationWorld_lib
  use Population_lib
  use mpi_lib
  implicit none

  type, public :: PopulationWorld

    type(Population) :: countries

  contains

    procedure :: give_and_take

  end type PopulationWorld

  interface PopulationWorld
    module procedure :: PopulationWorld_constructor
  end interface PopulationWorld

contains

  function PopulationWorld_constructor(n_generations, n_organisms, n_parents, n_genes, &
    n_tournament, n_crossover, fitness_fn) result(self)
    implicit none
    integer, intent(in) :: n_generations, n_organisms, n_parents, n_genes, n_tournament, n_crossover
    interface
      function fitness_fn(DNA) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: DNA
        real(kind=8) :: fit
      end function fitness_fn
    end interface
    type(PopulationWorld) :: self

    self%countries = Population(n_generations, n_organisms, n_parents, n_genes, &
    n_tournament, n_crossover, fitness_fn)

  end function PopulationWorld_constructor

  subroutine give_and_take(self)
    implicit none
    class(PopulationWorld), intent(in out) :: self

    integer :: i, status(MPI_STATUS_SIZE)
    type(Organism), dimension(10) :: immigrants, emigrants


    emigrants = self%countries%organisms(1:size(emigrants))

    do i=1, size(emigrants)

      if(rank.eq.0) then

        call MPI_Sendrecv(&
        emigrants(i)%DNA , emigrants%n_genes, MPI_DOUBLE_PRECISION, rank+1 , 0,
        immigrants(i)%DNA, emigrants%n_genes, MPI_DOUBLE_PRECISION, nproc-1, 1, &
        MPI_COMM_WORLD, status)

      else if(rank.eq.nproc-1) then

        call MPI_Sendrecv(&
        emigrants(i)%DNA , emigrants%n_genes, MPI_DOUBLE_PRECISION, 0     , 1,
        immigrants(i)%DNA, emigrants%n_genes, MPI_DOUBLE_PRECISION, rank-1, 0, &
        MPI_COMM_WORLD, status)

      else

        call MPI_Sendrecv(&
        emigrants(i)%DNA , emigrants%n_genes, MPI_DOUBLE_PRECISION, rank+1, 1,
        immigrants(i)%DNA, emigrants%n_genes, MPI_DOUBLE_PRECISION, rank-1, 0, &
        MPI_COMM_WORLD, status)

      end if

    end do

    self%countries%organisms(1:size(immigrants)) = imm

  end subroutine give_and_take


end module PopulationWorld_lib
