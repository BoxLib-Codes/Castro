module eos_override_module

  implicit none

  public eos_override

contains

  ! This is a user hook to override the details of the EOS state.

  subroutine eos_override(state)

    use eos_type_module, only: eos_t
    use actual_eos_module, only: eos_name
    use network
    use actual_network

    implicit none

    type (eos_t) :: state

    !$gpu

#ifdef NUCLEAR_EOS
    call fill_nuclear_state(state)
#endif

  end subroutine eos_override

#ifdef NUCLEAR_EOS
  subroutine fill_nuclear_state(eos_state)

    use network
    use eos_type_module, only: eos_t
    use UnitsModule, only: Gram, AtomicMassUnit

    implicit none

    type (eos_t), intent(inout) :: eos_state

    !$gpu

    eos_state % y_e = eos_state % xn(iye)
    eos_state % xne = eos_state % y_e * Gram * eos_state % rho / AtomicMassUnit
    eos_state % aux(ine) = eos_state % xne
    
  end subroutine fill_nuclear_state
#endif

end module eos_override_module
