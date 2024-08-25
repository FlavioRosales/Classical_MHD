subroutine overall_diagnosis
  use global_numbers
  implicit none

  fourier_space_gas = .false.

  mgas = gas%mass(fourier_space_gas)
  Kgas = gas%K_energy(fourier_space_gas)
  Wgas = gas%W_energy(fourier_space_gas)
  Ugas = gas%U_energy(fourier_space_gas)
  Egas = Kgas + Wgas + Ugas
  Qgas = 2.0d0 * (Kgas + Ugas) + Wgas
  pxgas = gas%px_lineal()
  pygas = gas%py_lineal()
  pzgas = gas%pz_lineal()
  pgas = dsqrt((pxgas(1)+pxgas(2))**2 + (pygas(1)+pygas(2))**2 + (pzgas(1)+pzgas(2))**2)
  rhomaxgas = gas%rhomax()
  xcmgas = gas%mass_center() / mgas

end subroutine overall_diagnosis
