module probdata_module

      use amrex_fort_module, only : rt => amrex_real
      use network, only: nspec

      real(rt), save :: xn(nspec)
      real(rt), save :: centx, centy, centz

      real(rt), save :: ye_err
      real(rt), save :: ye_grad
      real(rt), save :: ye_grad_rel

      integer, save :: max_ye_err_lev
      integer, save :: max_ye_grad_lev
      integer, save :: max_ye_grad_rel_lev

end module probdata_module
