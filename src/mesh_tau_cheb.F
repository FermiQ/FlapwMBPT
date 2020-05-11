      subroutine mesh_tau_spl
      use manager_mod
      use models_mod
      use parallel_mod
      use units_mod
      implicit none
      integer :: i_tau
      real*8 :: xk,hx
      tau_mesh(0)=0.d0
      tau_mesh(n_tau)=betta_t
c ---------------------------------------------------------------------
      hx=2.d0/n_tau
      do i_tau=1,n_tau/2
        xk=hx*i_tau
        tau_mesh(i_tau)=xk**exp_tau_gw*betta_t/2.d0
        tau_mesh(n_tau-i_tau)=betta_t-tau_mesh(i_tau)
      enddo
      end