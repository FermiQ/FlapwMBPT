      subroutine uu_prepare_ac
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,k,i,k1,iatom,isort,ind,jbas,jnd,ibas,j,ind_k,kd,
     &           ndimb,ind0,k11,k110,gbs_number,k10,ka0,ka,ib,k0,
     &           nbas,nbasi,irl,ibnd
      real*8 :: v(3),v1(3),gtild(3)
      real*8, allocatable :: fifi_j(:,:,:,:,:,:)
      complex*16, allocatable :: tmp(:,:),s(:,:),tmp1(:,:),ff(:,:,:,:),
     &                           aw(:,:,:,:),zw(:,:,:,:)
      if(allocated(psps_ac)) deallocate(psps_ac)
      allocate(psps_ac(nbndf_bnd,nbndf,num_baryc_ac_max,
     &                 ndim3_k(me_k+1),nspin))
      allocate(zw(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(aw(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      v1=0.d0
      aw=(0.d0,0.d0)
      zw=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          ka=n3_mpi_k(me_k+1)+ind_k
          k=k_c_from_a(ka)
          if(k==0) cycle
          do ibnd=1,n_low_bnd(ka,ispin)
            ib=ind_bands_bnd(ibnd,ka,ispin)
            zw(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
            aw(:,ibnd,k,ispin)=ev_bnd(:,ib,ind_k,ispin)
          enddo
        enddo
        do k=npnt_c+1,nqdiv_c
          ka=k_a_from_c(k)
          k0=i_kref_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(zw(1,1,k,ispin),zw(1,1,k0,ispin),k_group_c(k),
     &               pnt_c(1,k),nfun,n_low_bnd(ka0,ispin),maxb,
     &               io_lem,lmb,indbasa,limlb,0)
          call sym_a(aw(1,1,k,ispin),aw(1,1,k0,ispin),ka,k_group_c(k),
     &               pnt_c(1,k),nbasmpw,n_low_bnd(ka0,ispin),nbask(ka0),
     &               indgb(1,ka),iplf_bk(1,ka0))
        enddo
      enddo
      if(nproc_k/=1) then
        call dgop(zw,2*nfun*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
        call dgop(aw,2*nbasmpw*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
      endif
      allocate(tmp(maxel,nbndf))
      allocate(tmp1(nbndtm,nbndf))
      allocate(s(nbndtm,nbndtm))
      allocate(ff(maxel,maxel,natom,num_kk_k_ac))
      psps_ac=(0.d0,0.d0)
      allocate(fifi_j(maxlfun,maxlfun,0:2*maxb,nspin,nspin,nsort))
      call fifi_j_prepare(fifi_j,v1,1)
      do ispin=1,nspin
        call wan_wan_pw_0(fifi_j,ff,v1,1,ispin,ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do i=1,num_baryc_ac(k)
            call zone1_number(pnt_baryc_ac(1,i,k),rb0_c,ndiv_c,k1)
            k1=index_k1_c(k1)    !  0 BZ
            k10=i_kref_c(k1)
            ka0=k_a_from_c(k10)
            kd=1
c ----------------- MT contribution ----------------------------------
            do iatom=1,natom
              isort=is(iatom)
              ndimb=lfunm(isort)
              ind0=io_lem(iatom)
              call zgemm('n','n',ndimb,n_bnd(k,ispin),ndimb,(1.d0,0.d0),
     &	               ff(1,1,iatom,kd),maxel,
     &                   z_bnd(ind0,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp,
     &                   maxel)
              call zgemm('c','n',n_low_bnd(ka0,ispin),n_bnd(k,ispin),
     &	               ndimb,(1.d0,0.d0),zw(ind0,1,k1,ispin),nfun,tmp,
     &                   maxel,(1.d0,0.d0),psps_ac(1,1,i,ind_k,ispin),
     &                   nbndf_bnd)
            enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
            gtild=pnt_baryc_ac(:,i,k)-pnt_c(:,k1)
c ------ We form the overlap matrix -----------------------------
            k11=k_a_from_c(k1)
            k110=i_kref(k11)
            nbas=nbask(k)/nrel
            nbasi=nbask(k110)/nrel
            do jbas=1,nbas
              jnd=indgb(jbas,k)
              do ibas=1,nbasi
                ind=indgb(ibas,k11)
                v=gbs(:,jnd)-gbs(:,ind)+gtild
                j=gbs_number(v)
                if(complex_ro) then
                  s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
                else
                  s(ibas,jbas)=dcmplx(sovr(j),0.d0)
                endif
              enddo
            enddo
            do irl=1,nrel
              call zgemm('n','n',nbasi,n_bnd(k,ispin),nbas,
     &	               (1.d0,0.d0),s,nbndtm,
     &                   ev_bnd(1+nbas*(irl-1),1,ind_k,ispin),
     &                   nbasmpw,(0.d0,0.d0),tmp1,nbndtm)
              call zgemm('c','n',n_low_bnd(ka0,ispin),n_bnd(k,ispin),
     &	               nbasi,(1.d0,0.d0),
     &                   aw(1+nbasi*(irl-1),1,k1,ispin),nbasmpw,
     &                   tmp1,nbndtm,(1.d0,0.d0),
     &                   psps_ac(1,1,i,ind_k,ispin),nbndf_bnd)
            enddo  !! over irl
          enddo   !! over i
        enddo   !! over ind_k
      enddo   !! over ispin
      deallocate(fifi_j,ff,tmp,s,tmp1,aw,zw)
      end
      
      
      
      subroutine uu_prepare_ac_old
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,k,i,k1,iatom,isort,ind,jbas,jnd,ibas,j,ind_k,kd,
     &           ndimb,ind0,k11,k110,gbs_number,k10,ka0,ka,ib,k0,
     &           nbas,nbasi,irl,ibnd
      real*8 :: v(3),gtild(3)
      real*8, allocatable :: fifi_j(:,:,:,:,:,:,:)
      complex*16, allocatable :: tmp(:,:),s(:,:),tmp1(:,:),ff(:,:,:,:),
     &                           aw(:,:,:,:),zw(:,:,:,:)
      if(allocated(psps_ac)) deallocate(psps_ac)
      allocate(psps_ac(nbndf_bnd,nbndf,num_baryc_ac_max,
     &                 ndim3_k(me_k+1),nspin))
      allocate(zw(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(aw(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      aw=(0.d0,0.d0)
      zw=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          ka=n3_mpi_k(me_k+1)+ind_k
          k=k_c_from_a(ka)
          if(k==0) cycle
          do ibnd=1,n_low_bnd(ka,ispin)
            ib=ind_bands_bnd(ibnd,ka,ispin)
            zw(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
            aw(:,ibnd,k,ispin)=ev_bnd(:,ib,ind_k,ispin)
          enddo
        enddo
        do k=npnt_c+1,nqdiv_c
          ka=k_a_from_c(k)
          k0=i_kref_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(zw(1,1,k,ispin),zw(1,1,k0,ispin),k_group_c(k),
     &               pnt_c(1,k),nfun,n_low_bnd(ka0,ispin),maxb,
     &               io_lem,lmb,indbasa,limlb,0)
          call sym_a(aw(1,1,k,ispin),aw(1,1,k0,ispin),ka,k_group_c(k),
     &               pnt_c(1,k),nbasmpw,n_low_bnd(ka0,ispin),nbask(ka0),
     &               indgb(1,ka),iplf_bk(1,ka0))
        enddo
      enddo
      if(nproc_k/=1) then
        call dgop(zw,2*nfun*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
        call dgop(aw,2*nbasmpw*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
      endif
      allocate(tmp(maxel,nbndf))
      allocate(tmp1(nbndtm,nbndf))
      allocate(s(nbndtm,nbndtm))
      allocate(ff(maxel,maxel,natom,num_kk_k_ac))
      psps_ac=(0.d0,0.d0)
      allocate(fifi_j(maxlfun,maxlfun,0:2*maxb,nspin,nspin,nsort,
     &                num_kk_k_ac))
      call fifi_j_prepare(fifi_j,dif_kk_k_ac,num_kk_k_ac)
      do ispin=1,nspin
        call wan_wan_pw_0(fifi_j,ff,dif_kk_k_ac,num_kk_k_ac,ispin,ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do i=1,num_baryc_ac(k)
            call zone1_number(pnt_baryc_ac(1,i,k),rb0_c,ndiv_c,k1)
            k1=index_k1_c(k1)    !  0 BZ
            k10=i_kref_c(k1)
            ka0=k_a_from_c(k10)
            kd=ind_kk_k_ac(i,k)
c ----------------- MT contribution ----------------------------------
            do iatom=1,natom
              isort=is(iatom)
              ndimb=lfunm(isort)
              ind0=io_lem(iatom)
              call zgemm('n','n',ndimb,n_bnd(k,ispin),ndimb,(1.d0,0.d0),
     &	               ff(1,1,iatom,kd),maxel,
     &                   z_bnd(ind0,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp,
     &                   maxel)
              call zgemm('c','n',n_low_bnd(ka0,ispin),n_bnd(k,ispin),
     &	               ndimb,(1.d0,0.d0),zw(ind0,1,k1,ispin),nfun,tmp,
     &                   maxel,(1.d0,0.d0),psps_ac(1,1,i,ind_k,ispin),
     &                   nbndf_bnd)
            enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
            gtild=pnt_baryc_ac(:,i,k)-pnt_c(:,k1)
c ------ We form the overlap matrix -----------------------------
            k11=k_a_from_c(k1)
            k110=i_kref(k11)
            nbas=nbask(k)/nrel
            nbasi=nbask(k110)/nrel
            do jbas=1,nbas
              jnd=indgb(jbas,k)
              do ibas=1,nbasi
                ind=indgb(ibas,k11)
                v=gbs(:,jnd)-gbs(:,ind)+gtild
                j=gbs_number(v)
                if(complex_ro) then
                  s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
                else
                  s(ibas,jbas)=dcmplx(sovr(j),0.d0)
                endif
              enddo
            enddo
            do irl=1,nrel
              call zgemm('n','n',nbasi,n_bnd(k,ispin),nbas,
     &	               (1.d0,0.d0),s,nbndtm,
     &                   ev_bnd(1+nbas*(irl-1),1,ind_k,ispin),
     &                   nbasmpw,(0.d0,0.d0),tmp1,nbndtm)
              call zgemm('c','n',n_low_bnd(ka0,ispin),n_bnd(k,ispin),
     &	               nbasi,(1.d0,0.d0),
     &                   aw(1+nbasi*(irl-1),1,k1,ispin),nbasmpw,
     &                   tmp1,nbndtm,(1.d0,0.d0),
     &                   psps_ac(1,1,i,ind_k,ispin),nbndf_bnd)
            enddo  !! over irl
          enddo   !! over i
        enddo   !! over ind_k
      enddo   !! over ispin
      deallocate(fifi_j,ff,tmp,s,tmp1,aw,zw)
      end