      subroutine vertex_sig_1_inf(sig_1_1,iter_s)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_s
      complex*16, intent(out) :: sig_1_1(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      call vertex_sig_1_nlc(sig_1_1,iter_s)
      call timel('**** VERTEX_SIG_1_INF finished *****')
      end