      subroutine vertex_p_vmt_dir(tmatr)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: tmatr(n_lem2_red,n_lem2_red,ndim3_nu,
     &                                ndimc_kk(me_k+1),nspin,nspin)
      integer :: k,ispin,iq,ind_nu,i_nu,i,s2,np,ind_q,kq,ind_tau,ind_s2,
     &           ka,j
      real*8 :: p_herm,check_hermit,ph,v(3),gtild(3,2)
      real*8, allocatable :: poli(:,:),g_tau(:,:,:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),fif_full(:,:,:,:,:)
      complex*16, allocatable :: k0_pw(:,:,:,:,:),pol(:,:,:,:),
     &                           k_rs(:,:,:),z_red(:,:,:,:),
     &                           ps3(:,:,:,:),vrt(:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           z_full(:,:,:,:),k0_nu(:,:,:,:),
     &                           a_full(:,:,:,:),phf(:,:)
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(ax,phf,a_red)
c ------------------------------------------------------------------
      call norma_wf_red
c ------ Here for zero-order vertex in full basis ------------------
      allocate(fif_full(nrel,maxel**2,n_pbmtm_red,nsort,nspin))
      allocate(z_full(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_full(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      call aux_full(fif_full,z_full,a_full)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
      deallocate(g_omega)
      allocate(pol(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &             ndim_k_red(me_k+1)))
      pol=(0.d0,0.d0)
      allocate(ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,nqdiv_c))
      do ispin=1,nspin
        do ind_q=1,ndim_k_red(me_k+1)
          iq=n_mpi_k_red(me_k+1)+ind_q
          ka=k_a_from_c(iq)
          np=n_pbmt_red+nplwgw_red(ka)
          ps3=(0.d0,0.d0)
          do k=1,nqdiv_c
            v=pnt_c(:,k)-pnt_c(:,iq)
            call zone1_number(v,rb0_c,ndiv_c,kq)
            kq=index_k1_c(kq)    !! for K-Q
            gtild(:,1)=0.d0
            gtild(:,2)=v-pnt_c(:,kq)
            call ppm_factors_full(k,kq,iq,fif_full(1,1,1,1,ispin),
     &                            z_full(1,1,1,ispin),
     &                            a_full(1,1,1,ispin),ps3(1,1,1,k),
     &                            gtild,n_low_bnd(1,ispin))
          enddo   !! over k
          allocate(k_rs(n_lem2_red,np,ndim3_nu))
          k_rs=(0.d0,0.d0)
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            allocate(k0_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,2,ndim3_tau))
            do ind_tau=1,ndim3_tau
              call k0_bnd_stat_tau(ispin,ind_tau,iq,
     &                             k0_pw(1,1,1,1,ind_tau),
     &                             g_tau(1,1,1,1,1,ispin),ps3,s2)
            enddo
            allocate(k0_nu(nbndf_bnd,nbndf_bnd,nqdiv_c,ndim3_nu))
            call tau_to_nu_gen_spl(k0_nu,k0_pw,nbndf_bnd**2*nqdiv_c)
            deallocate(k0_pw)
c ------------------------------------------------------------------
            do ind_nu=1,ndim3_nu
              call k_real_space_vmt(iq,k0_nu(1,1,1,ind_nu),
     &                              k_rs(1,s2,ind_nu),z_red,ispin)
            enddo
            deallocate(k0_nu)
          enddo  !! ind_s2
          if(nproc_b/=1) call dgop(k_rs,2*n_lem2_red*np*ndim3_nu,'  +',
     &                             comm_b)
          allocate(vrt(n_lem2_red))
          do ind_s2=1,ndim4_pbr(me_b+1,iq)
            s2=n4_mpi_pbr(me_b+1,iq)+ind_s2
            do ind_nu=1,ndim3_nu
              call zgemm('n','n',n_lem2_red,1,n_lem2_red,(1.d0,0.d0),
     &                   tmatr(1,1,ind_nu,ind_q,ispin,ispin),n_lem2_red,
     &                   k_rs(1,s2,ind_nu),n_lem2_red,(0.d0,0.d0),vrt,
     &                   n_lem2_red)
              do i=1,np
                do j=1,n_lem2_red
                  pol(i,s2,ind_nu,ind_q)=pol(i,s2,ind_nu,ind_q)
     &              -conjg(k_rs(j,i,ind_nu))*vrt(j)
                enddo
              enddo
            enddo   !! over ind_nu
          enddo   !! over ind_s2
          deallocate(vrt,k_rs)
        enddo  !! over ind_q
      enddo   !! over ispin
c ---------------------------------------------------------------
      if(nproc_b/=1) call dgop(pol,
     &                2*n_pbtot_red**2*ndim3_nu*ndim_k_red(me_k+1),
     &                         '  +',comm_b)
      if(nspin==1.and.irel/=2) pol=2.d0*pol
      p_herm=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        iq=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(iq)
        np=n_pbmt_red+nplwgw_red(ka)
        allocate(poli(np,np))
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          ph=check_hermit(pol(1,1,ind_nu,ind_q),np,n_pbtot_red)
          p_herm=max(p_herm,ph)
          call pack_hermit(pol(1,1,ind_nu,ind_q),poli,np,
     &                     n_pbtot_red,np,0.d0,1.d0)
          p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)=
     &      p_vertex_nu_q(1:np,1:np,ind_nu,ind_q)+poli
        enddo
        deallocate(poli)
      enddo
      if(nproc_t/=1) call dgop(p_herm,1,'max',comm_t)
      if(nproc_k/=1) call dgop(p_herm,1,'max',comm_k)
      if(maswrk) write(iun,*)' VERTEX_P_VMT_DIR : P_HERM = ',p_herm
      call output_p_bnd_nu(2,pol,n_pbtot_red)
      deallocate(ps3,z_red,nind_fif,ind_fif,fif0,pol)
      call timel('**** VERTEX_P_VMT_DIR finished *****')
      end
