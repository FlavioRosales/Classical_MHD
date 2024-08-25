module Population_lib
  use Organism_lib
  implicit none

  type, public :: Population

    integer :: n_generations, n_organisms, n_parents, n_childs, n_tournament, n_crossover
    type(Organism) :: best_organism
    type(Organism), allocatable, dimension(:) :: organisms

  contains

    procedure :: classify
    procedure :: differential_mutation
    procedure :: selection_tournament
    procedure :: next_generation
    procedure :: elitism
    procedure :: evolve

  end type Population

  interface Population
    module procedure :: Population_constructor
  end interface Population

contains

  function Population_constructor(n_generations, n_organisms, n_parents, n_genes, &
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
    type(Population) :: self

    integer :: i
    real(kind=8), dimension(n_organisms) :: mutation_rates

    self%n_generations = n_generations
    self%n_organisms = n_organisms
    self%n_parents = n_parents
    self%n_childs = n_organisms - n_parents
    self%n_tournament = n_tournament
    self%n_crossover = n_crossover
    allocate(self%organisms(n_organisms))

    do i=1, n_organisms
      self%organisms(i) = Organism(n_genes = n_genes, fitness_fn = fitness_fn)

      if(rank.eq.master) print*, i, self%organisms(i)%fitness
    end do

    call self%classify

  end function Population_constructor

  recursive subroutine classify(self)
    implicit none
    class(Population), intent(in out) :: self

    integer :: i, n_left, n_right, n_repeating
    type(Organism), allocatable, dimension(:) :: organisms_left, organisms_right, organisms_repeating
    type(Population) :: left, right, repeating

    if(self%n_organisms>1) then

      n_left = 0
      n_right = 0
      n_repeating = 0
      do i=1, self%n_organisms
        if(self%organisms(i)%fitness>self%organisms(1)%fitness) then
          n_left = n_left + 1
        else if(self%organisms(i)%fitness<self%organisms(1)%fitness) then
          n_right = n_right + 1
        else
          n_repeating = n_repeating + 1
        end if
      end do

      allocate(organisms_left(n_left), organisms_right(n_right), organisms_repeating(n_repeating))

      n_left = 0
      n_right = 0
      n_repeating = 0
      do i=1, self%n_organisms
        if(self%organisms(i)%fitness>self%organisms(1)%fitness) then
          n_left = n_left + 1
          organisms_left(n_left) = self%organisms(i)
        else if(self%organisms(i)%fitness<self%organisms(1)%fitness) then
          n_right = n_right + 1
          organisms_right(n_right) = self%organisms(i)
        else
          n_repeating = n_repeating + 1
          organisms_repeating(n_repeating) = self%organisms(i)
        end if
      end do

      left%n_organisms = n_left
      left%organisms = organisms_left
      if(n_left>1) call left%classify
      deallocate(organisms_left)

      right%n_organisms = n_right
      right%organisms = organisms_right
      if(n_right>1) call right%classify
      deallocate(organisms_right)

      repeating%n_organisms = n_repeating
      repeating%organisms = organisms_repeating
      deallocate(organisms_repeating)

      do i=1, n_left
        self%organisms(i) = left%organisms(i)
      end do

      do i=n_left+1, n_left+n_repeating
        self%organisms(i) = repeating%organisms(i-n_left)
      end do

      do i=n_left+n_repeating+1, n_left+n_repeating+n_right
        self%organisms(i) = right%organisms(i-n_left-n_repeating)
      end do

    end if

  end subroutine classify

  subroutine differential_mutation(self, fitness_fn)
    implicit none
    class(Population), intent(in out) :: self
    interface
      function fitness_fn(DNA) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: DNA
        real(kind=8) :: fit
      end function fitness_fn
    end interface

    integer :: i
    integer, dimension(2) :: index
    real(kind=8), dimension(2) :: index_dp
    type(Organism) :: individuals(2), individual

    do i=1, self%n_organisms
      call random_mpi(index_dp)
      index = 1 + INT(index_dp * self%n_organisms)
      individuals(1) = self%organisms(index(1))
      individuals(2) = self%organisms(index(2))

      individual = self%organisms(i)
      individual%DNA = individual%DNA + 0.10d0 * (individuals(2)%DNA - individuals(1)%DNA)
      individual%fitness = fitness_fn(DNA = individual%DNA)

      if(individual%fitness>self%organisms(i)%fitness) self%organisms(i) = individual

    end do

  end subroutine differential_mutation

  subroutine selection_tournament(self, fitness_fn)
    implicit none
    class(Population), intent(in out) :: self
    interface
      function fitness_fn(x) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8) :: fit
      end function fitness_fn
    end interface

    integer, dimension(self%n_tournament) :: i_rivals
    real(kind=8), dimension(self%n_tournament) :: random_values
    type(Organism), dimension(self%n_tournament) :: rivals
    integer :: i, j, i_winner
    type(Organism), dimension(self%n_parents) :: winners

    do i = 1, self%n_parents
      ! Inicializa el ganador con el primer rival
      i_winner = 1

      ! Selección aleatoria de rivales
      call random_mpi(random_values)
      i_rivals = 1 + int(random_values * (self%n_organisms))

      ! Inicializa a los rivales seleccionados
      rivals = self%organisms(i_rivals)

      ! Encuentra al ganador dentro del torneo
      do j = 2, self%n_tournament
        if (rivals(j)%fitness > rivals(i_winner)%fitness) then
          i_winner = j
        end if
      end do

      ! Asigna al ganador como uno de los padres para la reproducción
      winners(i) = rivals(i_winner)

    end do

    self%organisms(1:self%n_parents) = winners

  end subroutine selection_tournament

  subroutine next_generation(self, fitness_fn)
    implicit none
    class(Population), intent(in out) :: self
    interface
      function fitness_fn(x) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8) :: fit
      end function fitness_fn
    end interface

    integer :: i
    integer, dimension(self%n_crossover) :: index
    real(kind=8), dimension(self%n_crossover) :: index_dp

    do i=1, self%n_childs

      call random_mpi(index_dp)
      index = 1 + INT(index_dp * self%n_parents)
      self%organisms(self%n_parents+i) = crossover(parents = self%organisms(index), fitness_fn = fitness_fn)

    end do

  end subroutine next_generation

  subroutine elitism(self, fitness_fn)
    implicit none
    class(Population), intent(in out) :: self
    interface
      function fitness_fn(x) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8) :: fit
      end function fitness_fn
    end interface

    call self%classify
    call self%next_generation(fitness_fn = fitness_fn)

  end subroutine elitism

  subroutine evolve(self, fitness_fn, filename)
    implicit none
    class(Population), intent(in out) :: self
    interface
      function fitness_fn(x) result(fit)
        implicit none
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8) :: fit
      end function fitness_fn
    end interface
    character(len=*), intent(in) :: filename

    integer :: i, j
    real(kind=8) :: error, fitness

    open(123, file = filename)

      do j=1, self%n_generations

        !call self%selection_tournament(fitness_fn)
        !call self%next_generation(fitness_fn)
        call self%elitism(fitness_fn = fitness_fn)
        call self%differential_mutation(fitness_fn)

        fitness = maxval(self%organisms(:)%fitness)
        error = 1.0d0 / fitness

        do i=1, self%n_organisms
          if(fitness.eq.self%organisms(i)%fitness) self%best_organism = self%organisms(i)
        end do

        if(rank.eq.master) write(*,*) 'Generation:', j, ' error:', error
        if(rank.eq.master) write(123,*) j, error, self%best_organism%DNA
      end do

    close(123)

  end subroutine evolve

end module Population_lib
