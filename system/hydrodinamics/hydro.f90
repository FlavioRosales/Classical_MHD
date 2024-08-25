module hydro_lib
  use hydro_initial_data_lib
  implicit none

  type, extends(hydro_initial_data) :: hydro
  contains

    final :: destroy

  end type hydro

  interface hydro
    module procedure :: hydro_constructor
  end interface hydro

contains

  !============================================================================!
  function hydro_constructor(K_poly, gamma, floor, EoS, flux_formula, boundary) result(this)
    implicit none
    real(kind=8), intent(in) :: K_poly, gamma, floor
    character(len=*), intent(in) :: EoS, flux_formula, boundary
    type(hydro) :: this

    if(rank.eq.master) print*, 'Creating a hydro variable...'

    if(.not.(EoS.eq.'ideal gas' .or. EoS.eq.'poly' .or. EoS.eq.'dust')) then
      stop "EoS can only be 'ideal gas', 'poly' or 'dust'."
    end if

    if(.not.(flux_formula.eq.'HLLE' .or. flux_formula.eq.'MARQUINA')) then
      stop "The flux formula can only be 'HLLE' or 'MARQUINA'."
    end if

    allocate(character(len = len(EoS)) :: this%EoS)
    this%EoS = EoS;

    allocate(character(len = len(flux_formula)) :: this%flux_formula)
    this%flux_formula = flux_formula

    this%gamma = gamma; this%floor = floor; this%K_poly = K_poly

    call this%hydro_memory
    allocate(character(len = len(boundary))  :: this%boundary )
    this%boundary = boundary

  end function hydro_constructor
  !============================================================================!
  subroutine destroy(this)
    implicit none
    type(hydro) :: this

    deallocate(this%rho, this%vx, this%vy, this%vz, this%e, this%p, this%V)
    deallocate(this%u, this%Fx, this%Fy, this%Fz, this%S)

  end subroutine destroy
  !============================================================================!

end module hydro_lib
