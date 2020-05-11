      subroutine k_t_sig(ispin,i_tau,k1,k_pw,z_red,a_red,ax,phf,grsmm,
     &                   grsmi,grsim,grsii,t_pw)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin,i_tau,k1
      real*8, intent(in) :: grsmm(nrel,nfun_red,nfun_red,0:n_tau,
     &                            nqdiv_c),
     &                      grsmi(nrel,nfun_red,nrel*nr_full_red,
     &                            0:n_tau,nqdiv_c),
     &                      grsim(nrel,nrel,nr_full_red,nfun_red,
     &                            0:n_tau,nqdiv_c),
     &                      grsii(nrel,nrel*nr_red_red,
     &                            nrel*nr_full_red,0:n_tau,nqdiv_c)
      complex*16, intent(in) :: k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      complex*16, allocatable :: t_rs(:,:,:),k_rs(:,:,:),krs0(:,:,:)
      integer :: i1_tau
      i1_tau=n_tau-i_tau
c ------------------------ MM ---------------------------------------
      allocate(k_rs(nfun_red,n_pbmt_red,nqdiv_c))
      call k_rs_sig_mm(k_pw,k_rs,z_red,ispin)
      allocate(t_rs(n_pbmt_red,nfun_red,nqdiv_c))
      if(irel/=2) then
        call tsig_mm(grsmm,k_rs,t_rs,fif0(1,1,1,1,1,ispin),i1_tau)
      else
        call tsig_mm_c(grsmm,k_rs,t_rs,fif0,i1_tau)
      endif
      deallocate(k_rs)
      call tsig_rr_to_k(n_pbmt_red,nfun_red,t_rs)
      call tsig_rs_to_ks_mm(t_rs,t_pw,k1,z_red,ispin)
      deallocate(t_rs)
c ----------------------- Common array for MI and II ----------------
      allocate(krs0(nbndf_bnd,nr_full_red,nqdiv_c))
      call k_rs_s_mi_ii(k1,k_pw,krs0,phf,ispin)
c ------------------------ MI ---------------------------------------
      allocate(k_rs(nfun_red,nr_full_red,nqdiv_c))
      call k_rs_s_mi(krs0,k_rs,z_red,ispin)
      allocate(t_rs(n_pbmt_red,nrel*nr_full_red,nqdiv_c))
      if(irel/=2) then
        call tsig_mi(grsmi,k_rs,t_rs,fif0(1,1,1,1,1,ispin),i1_tau)
      else
        call tsig_mi_c(grsmi,k_rs,t_rs,fif0,i1_tau)
      endif
      deallocate(k_rs)
      call tsig_rr_to_k(n_pbmt_red,nrel*nr_full_red,t_rs)
      call tsig_r_to_k_mi(t_rs,t_pw,k1,phf,ax,ispin)
      deallocate(t_rs)
c ------------------------ IM ---------------------------------------
      allocate(k_rs(nrel*nr_full_red,n_pbmt_red,nqdiv_c))
      call k_rs_s_im(k_pw,k_rs,a_red,phf,ispin)
      allocate(t_rs(nr_full_red,nfun_red,nqdiv_c))
      if(irel/=2) then
        call tsig_im(grsim,k_rs,t_rs,fif0(1,1,1,1,1,ispin),i1_tau)
      else
        call tsig_im_c(grsim,k_rs,t_rs,fif0,i1_tau)
      endif
      deallocate(k_rs)
      call tsig_rr_to_k(nr_full_red,nfun_red,t_rs)
      call tsig_r_to_k_im(t_rs,t_pw,k1,z_red,phf,ispin)
      deallocate(t_rs)
c ------------------------ II ---------------------------------------
      allocate(k_rs(nrel*nr_full_red,nr_full_red,nqdiv_c))
      call k_rs_s_ii(krs0,k_rs,a_red,phf,ispin)
      deallocate(krs0)
      allocate(t_rs(nr_full_red,nrel*nr_full_red,nqdiv_c))
      if(irel/=2) call tsig_ii(grsii,k_rs,t_rs,i1_tau)
      if(irel==2) call tsig_ii_c(grsii,k_rs,t_rs,i1_tau)
      deallocate(k_rs)
      call tsig_rr_to_k(nr_full_red,nrel*nr_full_red,t_rs)
      call tsig_r_to_k_ii(t_rs,t_pw,k1,ax,phf,ispin)
      deallocate(t_rs) 
      end
