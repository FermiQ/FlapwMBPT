      subroutine sig_1_stat_mt(it,ind_tau,t_pw,t_sing,g_tau,v2_q,ps3,s2,
     &                         k1,ispin,z_red,a_red,ax,phf,grsmm,
     &                         grsmi,grsim,grsii,sig_stat,nb,gctau)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: it,ind_tau,s2,k1,ispin,nb
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c),
     &                      grsmm(nrel,nfun_red,nfun_red,0:n_tau,
     &                            nqdiv_c),
     &                      grsmi(nrel,nfun_red,nrel*nr_full_red,
     &                            0:n_tau,nqdiv_c),
     &                      grsim(nrel,nrel,nr_full_red,nfun_red,
     &                            0:n_tau,nqdiv_c),
     &                      grsii(nrel,nrel*nr_red_red,
     &                            nrel*nr_full_red,0:n_tau,nqdiv_c)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c),
     &                          gctau(nbndf_bnd,nbndf_bnd,npnt_c,
     &                                 0:n_tau)
      complex*16, intent(inout) :: sig_stat(nbndf_bnd)
      complex*16, intent(out) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c),
     &                           t_sing(nbndf_bnd)
      integer :: i_tau,i1_tau,iq
      complex*16, allocatable :: k_pw(:,:,:),k_sing(:)
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      i1_tau=n_tau-i_tau
      allocate(k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c))
      call k_sig_stat(it,ind_tau,g_tau,ps3,v2_q,s2,k1,ispin,k_pw)
      allocate(k_sing(nbndf_bnd))
      call k_sig0_stat(it,g_tau(1,1,1,ind_tau,k1),s2,k1,ispin,k_sing)
      call k_t_sig(ispin,i_tau,k1,k_pw,z_red,a_red,ax,phf,grsmm,grsmi,
     &             grsim,grsii,t_pw)
      call k_t_sig_sing(ispin,k1,k_pw,k_sing,gctau(1,1,1,i1_tau),t_pw,
     &                  t_sing,ps3)
      deallocate(k_pw,k_sing)
      do iq=1,nqdiv_c
        call tsig_end_stat(k1,iq,t_pw(1,1,iq),g_tau,v2_q(1,1,iq),
     &                     ps3,sig_stat,ispin,nb,it,ind_tau)
      enddo  !! over iq
      call tsig0_end_stat(t_sing,gctau(1,1,k1,i_tau),sig_stat,nb)
      end
