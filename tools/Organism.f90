module Organism_lib
  use mpi_lib
  implicit none

  type, public :: Organism

    integer :: n_genes
    real(kind=8) :: fitness
    real(kind=8), allocatable, dimension(:) :: DNA

  contains

    procedure :: mutate

  end type Organism

  interface Organism
    module procedure :: Organism_constructor
  end interface Organism

  interface random_mpi
    module procedure :: random_mpi_const, random_mpi_vect
  end interface random_mpi

contains

  function Organism_constructor(n_genes, fitness_fn) result(self)
    implicit none
    integer, intent(in) :: n_genes
    interface
      function fitness_fn(DNA) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: DNA
        real(kind=8) :: fit
      end function fitness_fn
    end interface
    type(Organism) :: self

    self%n_genes = n_genes
    allocate(self%DNA(n_genes))
    call random_mpi(self%DNA)
    self%fitness = fitness_fn(DNA = self%DNA)

    print*, rank, self%fitness

    do while(self%fitness.ne.self%fitness)
      call random_mpi(self%DNA)
      self%fitness = fitness_fn(self%DNA)
    end do


  end function Organism_constructor

  subroutine mutate(self)
    implicit none
    class(Organism), intent(in out) :: self

    integer :: i, j
    real(kind=8) :: minprobmutate, mutation_rate(self%n_genes), mutation(self%n_genes)

    call random_mpi(minprobmutate)
    call random_mpi(mutation_rate)
    call random_mpi(mutation)

    minprobmutate = 0.750d0 + (1.0d0 - 0.75d0) * minprobmutate
    mutation_rate = -1.50d0 + (1.5d0 + 1.50d0) * mutation_rate

    do i=1, self%n_genes
      if(mutation(i) > minprobmutate) self%DNA(i) = mutation_rate(i) * self%DNA(i)
    end do

  end subroutine mutate

  function crossover(parents, fitness_fn) result(new_organism)
    implicit none
    type(Organism), dimension(:), intent(in) :: parents
    interface
      function fitness_fn(DNA) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: DNA
        real(kind=8) :: fit
      end function fitness_fn
    end interface
    type(Organism) :: new_organism

    integer :: i , j
    character(len=64) :: new_gen, gen, gen_parents(size(parents))
    integer, dimension(parents(1)%n_genes) :: id_parents
    real(kind=8), dimension(parents(1)%n_genes) :: random

    new_organism%n_genes = parents(1)%n_genes
    allocate(new_organism%DNA(new_organism%n_genes))

    call random_mpi(random)
    id_parents = 1 + int(random * size(parents))

    do i=1, new_organism%n_genes
      new_organism%DNA(i) = parents(id_parents(i))%DNA(i)
    end do

    call new_organism%mutate
    new_organism%fitness = fitness_fn(new_organism%DNA)

    do while(new_organism%fitness.ne.new_organism%fitness)
      call random_mpi(new_organism%DNA)
      new_organism%fitness = fitness_fn(new_organism%DNA)
    end do

  end function crossover

  character(len=64) function DoubleToBinary(inputNumber)
    implicit none
    real(kind=8), intent(in) :: inputNumber

    write(DoubleToBinary, '(b64.64)') inputNumber

  end function DoubleToBinary

  real(kind=8) function BinaryToDouble(inputBits)
    implicit none
    character(len=64), intent(in) :: inputBits

    read(inputBits, '(b64.64)') BinaryToDouble

  end function BinaryToDouble

  subroutine random_mpi_const(x)
    implicit none
    real(kind=8) :: x, x_local

    call random_number(x_local)
    call MPI_ALLREDUCE(x_local, x, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    x = x / dble(nproc)

  end subroutine random_mpi_const

  subroutine random_mpi_vect(x)
    implicit none
    real(kind=8), dimension(:) :: x
    real(kind=8), dimension(size(x)) :: x_local

    call random_number(x_local)

    call MPI_ALLREDUCE(x_local, x, size(x), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    x = x / dble(nproc)

  end subroutine random_mpi_vect

end module Organism_lib
