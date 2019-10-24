module probdata_module

      use amrex_fort_module, only : rt => amrex_real
      use network, only: nspec

      real(rt), save :: xn(nspec)
      real(rt), save :: rho_i, T_i, rhoe_i, p_i

end module probdata_module
