      subroutine vertex_sig_1_vmt(sig_1_1)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(out) :: sig_1_1(nbndf_bnd,nbndf_bnd,2,
     &                                   ndim3_tau,ndim_k_red(me_k+1),
     &                                   nspin)
      integer :: k,ispin,i,ii,ij,kl,l,mx2,i_tau,it,ind_tau,n,
     &           ind_k,ka,j,ir,i1_tau,ind_ir0,iatom,ir0,indi,indj,
     &           jsort,jatom,isort,n1,np,np1
      real*8 :: pi2,phase
      complex*16 :: cc1
      real*8, allocatable :: grsmm(:,:,:,:,:),w(:,:,:,:),w1(:,:,:,:),
     &                       aa(:,:,:,:),
     &                       g_tau(:,:,:,:,:,:),bb(:,:,:,:),
     &                       g_omega(:,:,:,:,:,:),cc(:,:,:,:),
     &                       b1(:,:,:,:),t(:,:),
     &                       vloc(:,:,:,:,:),sig_dyn(:,:,:,:,:,:)
      complex*16, allocatable :: z_red(:,:,:,:),sigt(:,:,:,:),
     &                           a_red(:,:,:,:),ax(:,:,:,:),
     &                           sig_stat(:,:,:,:,:),
     &                           sig_om(:,:),phf(:,:),sigr(:,:,:),
     &                           tmp(:,:)
      pi2=pi+pi
      mx2=maxel_red**2
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin))
      allocate(phf(nr_full_red,nqdiv_c))
      allocate(ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_bnd(z_red,a_red,phf,ax)
      deallocate(a_red,phf,ax,nind_fif,ind_fif,fif0)
c ---------- Interaction matrices in product basis -----------------
      allocate(vloc(mx2,mx2,natom,nspin,nspin))
      call w4_loc_stat(vloc)
c ------ Store G on OMEGA and TAU ---------------------------------
      allocate(g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,npnt_c,nspin))
      allocate(g_omega(nbndf_bnd,nbndf_bnd,2,ndim3_omega,npnt_c,nspin))
      call g_p_bnd_mm(g_tau,g_omega)
      deallocate(g_omega)
c -----------------
      allocate(grsmm(nrel,nfun_red,nfun_red,0:n_tau,nqdiv_c))
c ------------------------------------------------------------------
      allocate(sig_dyn(nfun_red,nfun_red,2,ndim3_tau,npnt_c,nspin))
      allocate(sig_stat(maxel_red,maxel_red,nrr_max_c,2,ndim3_tau))
      call timel('*** most allocations finished ******')
c ------------------------------------------------------------------
      sig_dyn=0.d0
      do ispin=1,nspin
        call g_mm_rs(g_tau(1,1,1,1,1,ispin),grsmm,z_red(1,1,1,ispin),
     &               ispin)
        call timel('**** G_REAL_SPACE finished *********')
        sig_stat=(0.d0,0.d0)
        do ind_ir0=1,ndim_nrr_c(me_k+1)
          ir0=n_mpi_nrr_c(me_k+1)+ind_ir0
          ir=nrr_red_c(1,ir0)
          iatom=nrr_red_c(2,ir0)
          jatom=nrr_red_c(3,ir0)
          jsort=is(jatom)
          isort=is(iatom)
          n=lfunm_red(isort)
          n1=lfunm_red(jsort)
          np=n_pbmt0_red(isort)
          np1=n_pbmt0_red(jsort)
          indi=io_lem_red(iatom)
          indj=io_lem_red(jatom)
          allocate(w(n,n,n,n))
          allocate(w1(n1,n1,n1,n1))
          kl=0
          do l=1,n
            do k=1,n
              kl=kl+1
              ij=0
              do j=1,n
                do i=1,n
                  ij=ij+1
                  w(i,j,k,l)=vloc(ij,kl,iatom,ispin,ispin)
                enddo
              enddo
            enddo
          enddo
          kl=0
          do l=1,n1
            do k=1,n1
              kl=kl+1
              ij=0
              do j=1,n1
                do i=1,n1
                  ij=ij+1
                  w1(i,j,k,l)=vloc(ij,kl,jatom,ispin,ispin)
                enddo
              enddo
            enddo
          enddo
c ------------------------------------------------------------------
          allocate(aa(n1,n,n,n))
          allocate(bb(n1,n1,n1,n))
          allocate(cc(n,n,n1,n1))
          allocate(b1(n,n1,n1,n1))
          allocate(t(n,n1))
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              do l=1,n
                do k=1,n
                  call dgemm('t','n',n1,n,n,1.d0,
     &                       grsmm(1,indi,indj,i_tau,ir),nfun_red,
     &                       w(1,1,k,l),n,0.d0,aa(1,1,k,l),n1)
                enddo
              enddo
              call dgemm('n','t',n1*n1*n1,n,n1,1.d0,w1,n1*n1*n1,
     &                   grsmm(1,indi,indj,i_tau,ir),nfun_red,0.d0,bb,
     &                   n1*n1*n1)
              do l=1,n1
                do k=1,n1
                  do j=1,n1
                    do i=1,n
                      b1(i,j,k,l)=bb(j,k,l,i)
                    enddo
                  enddo
                enddo
              enddo
              cc=0.d0
              do l=1,n1
                do k=1,n1
                  do j=1,n1
                    do i=1,n
                      do ii=1,n
                        cc(i,j,k,l)=cc(i,j,k,l)+aa(k,i,j,ii)
     &                          *grsmm(1,indi+ii-1,indj+l-1,i_tau,ir)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
              call dgemm('n','n',n,n1,n*n1*n1,1.d0,cc,n,b1,n*n1*n1,0.d0,
     &                   t,n)
              sig_stat(1:n,1:n1,ir0,it,ind_tau)=-t
            enddo  !!  over it
          enddo   !! over ind_tau
          deallocate(w,w1,aa,bb,cc,b1,t)
        enddo   !! over ind_ir0
        if(nproc_k/=1) then
          call dgop(sig_stat,4*maxel_red**2*ndim3_tau*nrr_max_c,
     &              '  +',comm_k)
        endif
        allocate(sigr(nfun_red,nfun_red,nqdiv_c))
        allocate(sigt(nfun_red,nfun_red,2,npnt))
        do ind_tau=1,ndim3_tau
          sigt=(0.d0,0.d0)
          do it=1,2
            if(irel/=2) call sym_sig_rs_c(sig_stat(1,1,1,it,ind_tau),
     &                                    sigr)
            if(irel==2) call sym_sig_rs_c_rel(sig_stat(1,1,1,it,
     &                                                 ind_tau),sigr)
            do k=1,npnt_c
              do ir=1,nqdiv_c
                phase=pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
                cc1=dcmplx(cos(phase),sin(phase))
                sigt(:,:,it,k)=sigt(:,:,it,k)+cc1*sigr(:,:,ir)
              enddo
            enddo
          enddo   !! over it
          do k=1,npnt_c
            call ferm_pack_tau(sigt(1,1,1,k),
     &                         sig_dyn(1,1,1,ind_tau,k,ispin),
     &                         nfun_red,nfun_red,nfun_red)
          enddo
        enddo
        deallocate(sigr,sigt)
        call timel('**** END of IR0 loop ***************')
      enddo   !! over ispin
      deallocate(g_tau,grsmm,vloc,sig_stat)
c ---------------------------------------------------------------
      allocate(sig_om(nfun_red,nfun_red))
      do ispin=1,nspin
        do ind_k=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_k
          ka=k_a_from_c(k)
          n=n_low_bnd(ka,ispin)
          allocate(tmp(n,nfun_red))
          do ind_tau=1,ndim3_tau
            do it=1,2
              call ferm_unpack_tau(sig_om,
     &                             sig_dyn(1,1,1,ind_tau,k,ispin),
     &                             nfun_red,nfun_red,nfun_red,it)
              call zgemm('c','n',n,nfun_red,nfun_red,(1.d0,0.d0),
     &                   z_red(1,1,k,ispin),nfun_red,sig_om,nfun_red,
     &                   (0.d0,0.d0),tmp,n)
              call zgemm('n','n',n,n,nfun_red,(1.d0,0.d0),tmp,n,
     &                   z_red(1,1,k,ispin),nfun_red,(0.d0,0.d0),
     &                   sig_1_1(1,1,it,ind_tau,ind_k,ispin),nbndf_bnd)
            enddo
          enddo
          deallocate(tmp)
        enddo
      enddo
      deallocate(sig_dyn,z_red,sig_om)
      call timel('**** VERTEX_SIG_1_VMT finished *****')
      end
