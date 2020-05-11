      subroutine add_sigma_vrt_bnd(ps)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ps
      integer :: ind_q,iq,k,ka,i,i0,j,j0,ind,ind_nu,k0,ka0,np,isy,ig,n,
     &           ind_k,ind_tau,ispin,kc,kst
      real*8 :: v(3),phase,pi2
      complex*16 :: bb,cc,cc1
      complex*16, allocatable :: p(:,:),dd(:),mm(:,:,:),mi(:,:,:),
     &                           ii(:,:,:),pmi(:,:),f(:),pii(:,:),
     &                           tmp(:,:,:),ab0(:,:,:),zb0(:,:,:),
     &                           ab(:,:),zb(:,:)
      e_vertex_sig_g_bnd=0.d0
      if(ps=='00000') return
      pi2=pi+pi
      allocate(dd(npnt))
      dd=(0.d0,0.d0)
      bb=(1.d0,0.d0)
      if(nqdiv_c/=1) then
        bb=bb+sum(bk_bz_red)
        do k=1,npnt
          dd(k)=sum(ckq_bz(:,k))
        enddo
      endif
      do ispin=1,nspin
        allocate(zb0(nfun,nbndf_bnd,npnt_c))
        allocate(ab0(nbasmpw,nbndf_bnd,npnt_c))
        zb0=(0.d0,0.d0)
        ab0=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          n=n_low_bnd(k,ispin)
          do i=1,n
            i0=ind_bands_bnd(i,k,ispin)
            zb0(:,i,kc)=z_bnd(:,i0,ind_k,ispin)
            ab0(:,i,kc)=ev_bnd(:,i0,ind_k,ispin)
          enddo
        enddo
        if(nproc_k/=1) then
          call dgop(zb0,2*nfun*nbndf_bnd*npnt_c,'  +',comm_k)
          call dgop(ab0,2*nbasmpw*nbndf_bnd*npnt_c,'  +',comm_k)
        endif
        do kc=1,npnt_c
          ka0=k_a_from_c(kc)
          n=n_low_bnd(ka0,ispin)
          allocate(tmp(n,n,2))
          do kst=1,k_star_c(kc)
            k=k_list_c(kst,kc)
            ka=k_a_from_c(k)
            allocate(zb(nfun,n))
            allocate(ab(nbasmpw,n))
            call sym_z(zb,zb0(1,1,kc),k_group_c(k),pnt_c(1,k),nfun,n,
     &                 maxb,io_lem,lmb,indbasa,limlb,0)
            call sym_a(ab,ab0(1,1,kc),ka,k_group_c(k),pnt_c(1,k),
     &                 nbasmpw,n,nbask(ka0),indgb(1,ka),
     &                 iplf_bk(1,ka0))
            do ind_tau=1,ndim3_tau
              call ferm_unpack_tau(tmp,sig_bnd_vertex(1,1,1,ind_tau,kc,
     &                                                ispin),n,n,
     &                             nbndf_bnd,3)
c ----------------------------- MM part ----------------------------
              do ind_q=1,ndim3_k(me_k+1)
                iq=n3_mpi_k(me_k+1)+ind_q
                if(k/=1) cc1=ckq_bz(k,iq)-bk_bz_red(k)
                if(k==1) cc1=bb-dd(iq)
                do j=1,n_pbmt_red
                  do i=1,n_pbmt_red
                    mm(i,j,ind_q)=mm(i,j,ind_q)+p(i,j)*cc1
                  enddo
                enddo
              enddo
c ----------------------------- MI part ----------------------------
c ----------------------------- II part ----------------------------
            enddo   !! over ind_tau
          enddo   !! over kst
        enddo   !! over kc
      enddo   !! over ispin
      deallocate(p,tmp,f,dd)
      call timel('*** ADD_SIGMA_VRT_BND finished *****')
      end
