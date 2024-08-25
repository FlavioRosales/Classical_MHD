module spherical_harmonics_lib
  use mesh_lib
  implicit none

  complex(kind=8), parameter :: i = (0.0d0, 1.0d0)
  real(kind=8), parameter :: pii = 1.0d0 / dacos(-1.0d0)

  function sph_harm(l,m) result(Y)
    implicit none
    complex(kind=8), dimension(domain%Nx, domain%Ny, domain%Nz) :: Y

    if(l.eq.0) then

      Y = 0.50d0 / dsqrt(pi)

    else if(l.eq.1) then

      if(m.eq.-1) then
        Y = + 0.50d0 * dsqrt(1.50d0 * pii) * exp(-i*domain%phi) * dsin(domain%theta)
      else if(m.eq.0) then
        Y = + 0.50d0 * dsqrt(3.0d0 * pii) * dcos(domain%theta)
      else if(m.eq.+1) then
        Y = - 0.50d0 * dsqrt(1.50d0 * pii) * exp(+i*domain%phi) * dsin(domain%theta)
      end if

    else if(l.eq.2) then

      if(m.eq.-2) then
        Y = + 0.250d0 * dsqrt(7.50d0 * pii) * exp(-2.0d0*i*domain%phi) * dsin(domain%theta)**2
      else if(m.eq.-1) then
        Y = + 0.250d0 * dsqrt(7.50d0 * pii) * exp(-1.0d0*i*domain%phi) * dsin(domain%theta) * dcos(domain%theta)
      else if(m.eq.0) then
        Y = + 0.250d0 * dsqrt(5.0d0 * pii) * (3.0d0 * dcos(domain%theta)**2 - 1.0d0)
      else if(m.eq.+1) then
        Y = - 0.250d0 * dsqrt(7.50d0 * pii) * exp(+1.0d0*i*domain%phi) * dsin(domain%theta) * dcos(domain%theta)
      else if(m.eq.-2) then
        Y = + 0.250d0 * dsqrt(7.50d0 * pii) * exp(-2.0d0*i*domain%phi) * dsin(domain%theta)**2
      end if

    else if(l.eq.3) then

      if(m.eq.-3) then
        Y = + 0.1250d0 * dsqrt(35.0d0 * pii) * exp(-3.0d0*i*domain%phi) * dsin(domain%theta)**3
      else if(m.eq.-2) then
        Y = + 0.2500d0 * dsqrt(52.5d0 * pii) * exp(-2.0d0*i*domain%phi) * dsin(domain%theta)**2 * dcos(domain%theta)
      else if(m.eq.-1) then
        Y = + 0.1250d0 * dsqrt(21.0d0 * pii) * exp(-1.0d0*i*domain%phi) * dsin(domain%theta) * (5.0d0 * dcos(domain%theta)**2 - 1.0d0)
      else if(m.eq.0) then
        Y = + 0.2500d0 * dsqrt(7.00d0 * pii) * (5.0d0 * dcos(domain%theta)**3 - 3.0d0 * dcos(domain%theta))
      else if(m.eq.+1) then
        Y = - 0.1250d0 * dsqrt(21.0d0 * pii) * exp(+1.0d0*i*domain%phi) * dsin(domain%theta) * (5.0d0 * dcos(domain%theta)**2 - 1.0d0)
      else if(m.eq.+2) then
        Y = + 0.2500d0 * dsqrt(52.5d0 * pii) * exp(+2.0d0*i*domain%phi) * dsin(domain%theta)**2 * dcos(domain%theta)
      else if(m.eq.+3) then
        Y = - 0.1250d0 * dsqrt(35.0d0 * pii) * exp(+3.0d0*i*domain%phi) * dsin(domain%theta)**3
      end if

    else if(l.eq.4) then

      if(m.eq.-4) then
        Y = + 0.18750d0 * dsqrt(17.5d0 * pii) * exp(-4.0d0*i*domain%phi) * dsin(domain%theta)**4
      else if(m.eq.-3) then
        Y = + 0.37500d0 * dsqrt(35.0d0 * pii) * exp(-3.0d0*i*domain%phi) * dsin(domain%theta)**3 * dcos(domain%theta)
      else if(m.eq.-2) then
        Y = + 0.37500d0 * dsqrt(2.500d0 * pii) * exp(-2.0d0*i*domain%phi) * dsin(domain%theta)**2 * (7.0d0 * dcos(domain%theta)**2 - 1.0d0)
      else if(m.eq.-1) then
        Y = + 0.37500d0 * dsqrt(5.000d0 * pii) * exp(-1.0d0*i*domain%phi) * dsin(domain%theta) * (7.0d0 * dcos(domain%theta)**3 - 3.0d0 * dcos(theta))
      else if(m.eq.0) then
        Y = + 0.18750d0 * dsqrt(1.000d0 * pii) * (35.0d0 * )
      else if(m.eq.+1) then
      else if(m.eq.+2) then
      else if(m.eq.+3) then
      else if(m.eq.+4) then
      end if

    end if

  end function sph_harm

end module spherical_harmonics_lib
