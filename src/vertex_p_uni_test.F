      subroutine vertex_p_uni_test(iter_w,iter_t)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iter_w,iter_t
      logical :: spin_flips
      integer :: k,ispin,ind_nu,i,it_vrt,nit,ib,ind_omega,i_tau,n,nk,
     &           kf,ka,j,kst,k0,iq
      real*8 :: de,green0,conv_1,pi2,pol
      complex*16 :: cc
      real*8, allocatable :: wnu(:,:,:),
     &                       w_tau(:,:,:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),ex(:,:,:),
     &                       g_x(:,:,:,:),
     &                       fif_full(:,:,:,:,:),
     &                       k0_pw(:,:,:,:,:,:),
     &                       k_pw(:,:,:,:,:,:),lambda_dyn(:,:,:,:,:,:),
     &                       vrt_dyn_old(:,:,:,:,:,:)
      complex*16, allocatable :: lambda_stat(:,:,:,:),aa(:,:),
     &                           ps3(:,:,:),v2_q(:,:,:),
     &                           vrt_stat_old(:,:,:,:),bb(:,:),
     &                           z_full(:,:,:,:),a_full(:,:,:,:)
      pi2=pi+pi
c ------ Here for zero-order vertex in full basis ------------------
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
      allocate(v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c))
      v2_q=(0.d0,0.d0)
      do k=1,npnt_c
        ka=k_a_from_c(k)
        do kst=1,k_star_c(k)
          iq=k_list_c(kst,k)
          v2_q(:,:,iq)=v_red_q(:,:,ka)
          call sym_w_red(iq,k,v2_q(1,1,iq))
        enddo   !! over kst
      enddo   !! over k
c ------ Store G on OMEGA and TAU ----------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
c ---------- Interaction W on TAU ----------------------------------
      allocate(w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,nqdiv_c))
      allocate(wnu(n_pbtot_red,n_pbtot_red,ndim3_nu))
      allocate(bb(n_pbtot_red,n_pbtot_red))
      do k=1,nqdiv_c
        k0=i_kref_c(k)
        ka=k_a_from_c(k0)
        n=n_pbmt_red+nplwgw_red(ka)
        do ind_nu=1,ndim3_nu
          bb=w_red_q(:,:,ind_nu,ka)
          call sym_w_red(k,k0,bb)
          call pack_hermit(bb,wnu(1,1,ind_nu),n,
     &                     n_pbtot_red,n_pbtot_red,0.d0,1.d0)
        enddo
        call nu_to_tau_ba(wnu,n_pbtot_red,w_tau(1,1,1,k),n_pbtot_red,n)
      enddo
      deallocate(wnu,bb)
c ---- Exchange eigen values relative to Chem_POT ------------------
      allocate(ex(nbndf_bnd,npnt_c,nspin))
      ex=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          kf=k_a_from_c(k)
          do i=1,n_low_bnd(kf,ispin)
            ib=ind_bands_bnd(i,kf,ispin)
            ex(i,k,ispin)=e_bnd(ib,kf,ispin)-chem_pot
          enddo
        enddo   !! over k
      enddo  !! over ispin
c --------- Exchange Green's function ------------------------------
      allocate(g_x(nbndf_bnd,0:n_tau,npnt_c,nspin))
      g_x=0.d0
      do ispin=1,nspin
        do k=1,npnt_c
          kf=k_a_from_c(k)
          do i=1,n_low_bnd(kf,ispin)
            de=ex(i,k,ispin)
            do i_tau=0,n_tau
              g_x(i,i_tau,k,ispin)=green0(de,tau_mesh(i_tau))
            enddo
          enddo
        enddo   !! over k
      enddo  !! over ispin
      allocate(k0_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(k_pw(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(lambda_stat(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
      allocate(vrt_stat_old(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
c ------------------------------------------------------------------
      pol=0.d0
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red))
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------
      call timel('**** VERTEX_BND_GW : loop starts ***')
      nit=iter_w
      allocate(lambda_dyn(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,
     &                    nspin))
      allocate(vrt_dyn_old(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,
     &                     nspin))
c ------------------------------------------------------------------
      call timel('**** K0_BND starts *****************')
      lambda_dyn=0.d0
c -------- Get zero order generalized susceptibility ---------------
      call k0_bnd_uni_test(k0_pw,g_omega,g_x,ex)
      call timel('**** K0_BND finished ***************')
c ------------------------------------------------------------------
      k_pw=0.d0
      vrt_dyn_old=0.d0
      vrt_stat_old=(0.d0,0.d0)
      do it_vrt=1,nit
        call timel('****** Iter-n starts ***************')
        k_pw=k0_pw+k_pw
        spin_flips=.false.
        if(it_vrt<=iter_t) spin_flips=.true.
        call vertex_wk_uni_test(fif_full,k_pw,g_tau,z_full,a_full,w_tau,
     &                          lambda_dyn,lambda_stat,v2_q)
c --------------------------------------------------------------------
        conv_1=0.d0
        do ispin=1,nspin
          do k=1,npnt_c
            ka=k_a_from_c(k)
            nk=n_low_bnd(ka,ispin)
            allocate(aa(nk,nk))
            allocate(bb(nk,nk))
            do ind_omega=1,ndim3_omega
              call ferm_unpack_omega(aa,
     &                    lambda_dyn(1,1,1,ind_omega,k,ispin),nk,
     &                               nbndf_bnd,nk)
              call ferm_unpack_omega(bb,
     &                    vrt_dyn_old(1,1,1,ind_omega,k,ispin),nk,
     &                               nbndf_bnd,nk)
              do j=1,nk
                do i=1,nk
                  cc=aa(i,j)+lambda_stat(i,j,k,ispin)
     &              -bb(i,j)-vrt_stat_old(i,j,k,ispin)
                  conv_1=max(conv_1,abs(cc))
                enddo
              enddo
            enddo
            deallocate(aa,bb)
          enddo
        enddo
        if(maswrk) write(iun,*)it_vrt,conv_1,' :Conv'
        vrt_dyn_old=lambda_dyn
        vrt_stat_old=lambda_stat
        call k_bnd_uni_test(k_pw,g_omega,lambda_dyn,lambda_stat,g_x,ex)
        if(it_vrt==nit.and.me_t==0) then
          do ispin=1,nspin
            do k=1,npnt_c
              ka=k_a_from_c(k)
              nk=n_low_bnd(ka,ispin)
              do i=1,nk
                cc=0.5d0*(k_pw(i,i,1,1,k,ispin)+k_pw(i,i,2,1,k,ispin))
                pol=pol-wgt_c(k)*cc/sqrt(amega)
              enddo
            enddo
          enddo
        endif
      enddo  !! over it_vrt
      call output_vrt_uni(lambda_stat,lambda_dyn)
      deallocate(lambda_dyn,vrt_dyn_old)
      deallocate(k0_pw,k_pw,lambda_stat,g_x,ex,vrt_stat_old,fif_full,
     &           z_full,a_full,g_tau,g_omega,w_tau,v2_q)
c ---------------------------------------------------------------
      if(nproc_t/=1) call dgop(pol,1,'  +',comm_t)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      if(maswrk) write(iun,*)' Uniform P:',pol
c ------------------------------------------------------------------
      call timel('**** VERTEX_P_UNI_TEST finished ****')
      end
