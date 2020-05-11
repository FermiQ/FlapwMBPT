      subroutine vertex_sig_inf(psi)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: psi
      integer :: ispin,ind_tau,ind_k,k,ka,n,it,mx2
      real*8, allocatable :: g_rr_tau(:,:,:,:,:,:),wst4loc(:,:,:,:,:)
      complex*16, allocatable :: sig_1_1(:,:,:,:,:,:),sig(:,:,:),
     &                           sig_1_2(:,:,:,:,:,:),
     &                           sig_4_2(:,:,:,:,:,:),
     &                           sig_56_2(:,:,:,:,:,:),
     &                           sig_psi5_ph(:,:,:,:,:,:),
     &                           sig_psi5_pp(:,:,:,:,:,:),
     &                           s(:,:,:,:,:,:)
      mx2=maxel_red**2
      if(psi=='00000') return
      if(allocated(sig_bnd_vertex)) deallocate(sig_bnd_vertex)
      allocate(sig_bnd_vertex(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,
     &                        nspin))
      sig_bnd_vertex=0.d0
      allocate(w_red_q(n_pbtot_red,n_pbtot_red,ndim3_nu,npnt))
      call w_c_from_b
c -------------------------------------------------------------------      
      if(psi(1:1)/='0') then
        allocate(sig_1_1(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                   ndim_k_red(me_k+1),nspin))
        call vertex_sig_1_inf(sig_1_1,iter_sigma_gwg)
      endif
      if(psi(2:2)/='0') then
        allocate(sig_56_2(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                    ndim_k_red(me_k+1),nspin))
        if(psi(2:2)=='5') call vertex_sig_2_56(sig_56_2,'sta')
        if(psi(2:2)=='6') call vertex_sig_2_56(sig_56_2,'dyn')
      endif
      if(psi(3:3)/='0') then
        allocate(sig_1_2(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                   ndim_k_red(me_k+1),nspin))
        allocate(sig_4_2(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                   ndim_k_red(me_k+1),nspin)) 
        call vertex_sig_2_14(sig_1_2)
      endif
      if(psi(5:5)/='0') then
        allocate(sig_psi5_ph(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                       ndim_k_red(me_k+1),nspin))
        allocate(sig_psi5_pp(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                       ndim_k_red(me_k+1),nspin))
        if(psi(5:5)=='3') then
          allocate(g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,ndim3_tau,
     &                      nspin))
          do ispin=1,nspin
            call g_real_space_mm(g_rr_tau(1,1,1,1,1,ispin),ispin)
          enddo
          allocate(wst4loc(mx2,mx2,natom,nspin,nspin))
          call w4_loc_stat(wst4loc)
          call sig_psi5_g2_w1_stat(sig_psi5_ph,g_rr_tau,wst4loc,'ph')
          call sig_psi5_g2_w1_stat(sig_psi5_pp,g_rr_tau,wst4loc,'pp')
          deallocate(g_rr_tau,wst4loc)
        endif
      endif
      allocate(sig(nbndf_bnd,nbndf_bnd,2))
      do ispin=1,nspin
        do ind_tau=1,ndim3_tau
          do ind_k=1,ndim_k_red(me_k+1)
            k=n_mpi_k_red(me_k+1)+ind_k
            ka=k_a_from_c(k)
            n=n_low_bnd(ka,ispin)
            if(psi(3:3)/='0') then
              do it=1,2
                sig_4_2(1:n,1:n,it,ind_tau,ind_k,ispin)=
     &            conjg(transpose(sig_1_2(1:n,1:n,it,ind_tau,ind_k,
     &                                    ispin)))
              enddo
            endif
            sig=(0.d0,0.d0)  
            do it=1,2
              if(psi(1:1)/='0') sig(1:n,1:n,it)=sig(1:n,1:n,it)
     &                     +sig_1_1(1:n,1:n,it,ind_tau,ind_k,ispin)
              if(psi(2:2)/='0') sig(1:n,1:n,it)=sig(1:n,1:n,it)
     &                       +sig_56_2(1:n,1:n,it,ind_tau,ind_k,ispin)
              if(psi(3:3)/='0') sig(1:n,1:n,it)=sig(1:n,1:n,it)
     &                       +sig_1_2(1:n,1:n,it,ind_tau,ind_k,ispin)
     &                       +sig_4_2(1:n,1:n,it,ind_tau,ind_k,ispin)
              if(psi(5:5)/='0') sig(1:n,1:n,it)=sig(1:n,1:n,it)
     &                    +sig_psi5_ph(1:n,1:n,it,ind_tau,ind_k,ispin)
     &                    +sig_psi5_pp(1:n,1:n,it,ind_tau,ind_k,ispin)
            enddo
            call ferm_pack_tau(sig,
     &                         sig_bnd_vertex(1,1,1,ind_tau,k,ispin),
     &                         n,nbndf_bnd,nbndf_bnd)
          enddo
        enddo
      enddo
      deallocate(sig)
      if(nproc_k/=1) call dgop(sig_bnd_vertex,
     &                         2*nbndf_bnd**2*ndim3_tau*npnt_c*nspin,
     &                         '  +',comm_k)
      if(psi(1:1)/='0') call output_sigma_vrt(sig_1_1,'FLL1',
     &                                            evolt/2)
      if(psi(2:2)/='0') call output_sigma_vrt(sig_56_2,'56_2',
     &                                            evolt/2)
      if(psi(3:3)/='0') then
        call output_sigma_vrt(sig_1_2,'_1_2',evolt/2)
        call output_sigma_vrt(sig_4_2,'_4_2',evolt/2)
      endif
      allocate(s(nbndf_bnd,nbndf_bnd,2,ndim3_tau,ndim_k_red(me_k+1),
     &           nspin))
      if(psi(5:5)/='0') then
        call output_sigma_vrt(sig_psi5_ph,'5_ph',evolt/2)
        call output_sigma_vrt(sig_psi5_pp,'5_pp',evolt/2)
        s=sig_psi5_ph+sig_psi5_pp
        call output_sigma_vrt(s,'FLL5',evolt/2)
      endif
      s=(0.d0,0.d0)  
      if(psi(2:4)/='000') then
        if(psi(2:2)/='0') s=s+sig_56_2  
        if(psi(3:3)/='0') s=s+sig_1_2+sig_4_2
        call output_sigma_vrt(s,'FLL2',evolt/2)
      endif
      if(psi(1:1)/='0') s=s+sig_1_1
      if(psi(5:5)/='0') then
        s=s+sig_psi5_ph+sig_psi5_pp
        deallocate(sig_psi5_ph,sig_psi5_pp)
      endif
      call output_sigma_vrt(s,'FULL',evolt/2)
      if(psi(1:1)/='0') deallocate(sig_1_1)
      if(psi(2:2)/='0') deallocate(sig_56_2)
      if(psi(3:3)/='0') deallocate(sig_1_2,sig_4_2)
      deallocate(s,w_red_q)
      end
