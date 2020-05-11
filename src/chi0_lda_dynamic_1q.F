      subroutine chi0_lda_dynamic_1q(iq,isp1,isp2,chi0,n,wrx,nd,ns)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: isp1,isp2,n,iq,nd,ns
      real*8, intent(in) :: wrx(0:nrax_chi)
      complex*16, intent(out) :: chi0(n,n,nd)
      integer :: k,ind_k,i,j,jj,ig,n1,n2,n12,kq,kq0,ii,k0,ks,ind_nu,
     &           i_nu,nq0
      real*8 :: gtild(3,2),de,v(3),fermi_dirac,fk,fkq,d,w,ek,ekq
      complex*16 :: nuc
      integer, allocatable :: indk0(:),ndk(:),nsk(:)
      real*8, allocatable :: gb(:,:),fifim(:,:,:,:),wk0(:),gre(:,:,:)
      complex*16, allocatable :: ps3(:,:,:),tmp(:,:,:),tmp1(:,:),
     &                           t(:,:,:),zb02(:,:,:),ev02(:,:,:),
     &                           zb2(:,:),ev2(:,:),zb01(:,:,:),
     &                           ev01(:,:,:),zb1(:,:),ev1(:,:),
     &                           aa(:,:,:),ea(:,:,:)
      allocate(tmp(n,n,nd))
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
      call fi0_full(fifim,isp1,isp2)
      allocate(gb(nbndf,2))
      chi0=(0.d0,0.d0)
      gtild(:,1)=0.d0
      allocate(zb02(nfun,nbndf,npnt))
      allocate(ev02(nbasmpw,nbndf,npnt))
      zb02=(0.d0,0.d0)
      ev02=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        zb02(:,:,k)=z_bnd(:,:,ind_k,isp2)
        ev02(:,:,k)=ev_bnd(:,:,ind_k,isp2)
      enddo
      if(isp1/=isp2) then
        allocate(zb01(nfun,nbndf,npnt))
        allocate(ev01(nbasmpw,nbndf,npnt))
        zb01=(0.d0,0.d0)
        ev01=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          zb01(:,:,k)=z_bnd(:,:,ind_k,isp1)
          ev01(:,:,k)=ev_bnd(:,:,ind_k,isp1)
        enddo
      endif
      if(nproc_k/=1) then
        call dgop(zb02,2*nfun*nbndf*npnt,'  +',comm_k)
        call dgop(ev02,2*nbasmpw*nbndf*npnt,'  +',comm_k)
        if(isp1/=isp2) then
          call dgop(zb01,2*nfun*nbndf*npnt,'  +',comm_k)
          call dgop(ev01,2*nbasmpw*nbndf*npnt,'  +',comm_k)
        endif
      endif
      allocate(zb2(nfun,nbndf))
      allocate(ev2(nbasmpw,nbndf))
      allocate(zb1(nfun,nbndf))
      allocate(ev1(nbasmpw,nbndf))
c ---------------------------------------------------------------
      allocate(indk0(nqdiv))
      allocate(wk0(nqdiv))
      call small_star_q(iq,nq0,indk0,wk0)
      allocate(ndk(nproc_k))
      allocate(nsk(nproc_k))
      call size_shift_par(nq0,nproc_k,ndk,nsk)
      tmp=(0.d0,0.d0)
      do ind_k=1,ndk(me_k+1)
        ks=nsk(me_k+1)+ind_k
        k=indk0(ks)
        k0=i_kref(k)
        v=pnt(:,k)-pnt(:,iq)
        call zone1_number(v,rb0,ndiv,kq)
        kq=index_k1(kq)    !! for K-Q
        kq0=i_kref(kq)
        gtild(:,2)=v-pnt(:,kq)
        n1=n_bnd(k0,isp1)
        n2=n_bnd(kq0,isp2)
        do i=1,n1
          de=e_bnd(i,k0,isp1)-chem_pot
          gb(i,1)=fermi_dirac(de)
        enddo
        do i=1,n2
          de=e_bnd(i,kq0,isp2)-chem_pot
          gb(i,2)=fermi_dirac(de)
        enddo
        allocate(ps3(n1,n2,n))
        call sym_z(zb2,zb02(1,1,kq0),k_group(kq),pnt(1,kq),nfun,n2,maxb,
     &             io_lem,lmb,indbasa,limlb,0)
        call sym_a(ev2,ev02(1,1,kq0),kq,k_group(kq),pnt(1,kq),nbasmpw,
     &             n2,nbask(kq0),indgb(1,kq),iplf_bk(1,kq0))
        if(isp1==isp2) then
          call sym_z(zb1,zb02(1,1,k0),k_group(k),pnt(1,k),nfun,n2,maxb,
     &               io_lem,lmb,indbasa,limlb,0)
          call sym_a(ev1,ev02(1,1,k0),k,k_group(k),pnt(1,k),nbasmpw,
     &               n2,nbask(k0),indgb(1,k),iplf_bk(1,k0))
        else
          call sym_z(zb1,zb01(1,1,k0),k_group(k),pnt(1,k),nfun,n2,maxb,
     &               io_lem,lmb,indbasa,limlb,0)
          call sym_a(ev1,ev01(1,1,k0),k,k_group(k),pnt(1,k),nbasmpw,
     &               n2,nbask(k0),indgb(1,k),iplf_bk(1,k0))
        endif
        call ppm_factors0(k,kq,iq,fifim,zb1,zb2,ev1,ev2,ps3,gtild,n,n1,
     &                    n2,0,nbndf)
        allocate(t(n1,n2,n))
        allocate(aa(n1,n2,nd))
        if(nn_k_int==1) then
          do j=1,n2
            do i=1,n1
              de=e_bnd(j,kq0,isp2)-e_bnd(i,k0,isp1)
              d=gb(j,2)-gb(i,1)
              do ind_nu=1,nd
                i_nu=ns+ind_nu-1
                nuc=dcmplx(wrx(i_nu),e_small_bos)
                aa(i,j,ind_nu)=d/(nuc+de)
              enddo
            enddo
          enddo
          aa=aa*wk0(ks)
        else
          w=wk0(ks)/nn_k_int
          allocate(gre(nbndf,3,2))
          allocate(ea(n1,n1,3))
          call bands_gradient(nbndf,n1,k,isp1,zb1,ev1,ea)
          do i=1,n1
            gre(i,:,1)=(0.d0,-2.d0)*ea(i,i,:)
          enddo
          deallocate(ea)
          allocate(ea(n2,n2,3))
          call bands_gradient(nbndf,n2,kq,isp2,zb2,ev2,ea)
          do i=1,n2
            gre(i,:,2)=(0.d0,-2.d0)*ea(i,i,:)
          enddo
          deallocate(ea)
          aa=(0.d0,0.d0)
          do j=1,n2
            do i=1,n1
              do ii=1,nn_k_int
                ekq=e_bnd(j,kq0,isp2)+gre(j,1,2)*k_int_vec(1,ii)
     &                               +gre(j,2,2)*k_int_vec(2,ii)
     &                               +gre(j,3,2)*k_int_vec(3,ii)
                ek=e_bnd(i,k0,isp1)+gre(j,1,1)*k_int_vec(1,ii)
     &                             +gre(j,2,1)*k_int_vec(2,ii)
     &                             +gre(j,3,1)*k_int_vec(3,ii)
                de=ekq-ek
                fkq=fermi_dirac(ekq-chem_pot)
                fk=fermi_dirac(ek-chem_pot)
                d=w*(fkq-fk)
                do ind_nu=1,nd
                  i_nu=ns+ind_nu-1
                  nuc=dcmplx(wrx(i_nu),e_small_bos)
                  aa(i,j,ind_nu)=aa(i,j,ind_nu)+d/(nuc+de)
                enddo
              enddo
            enddo
          enddo
          deallocate(gre)
        endif
        n12=n1*n2
        do ind_nu=1,nd
          do j=1,n2
            do i=1,n1
              do jj=1,n
                t(i,j,jj)=aa(i,j,ind_nu)*ps3(i,j,jj)
              enddo
            enddo
          enddo
          call zgemm('c','n',n,n,n12,(1.d0,0.d0),ps3,n12,t,n12,
     &               (1.d0,0.d0),tmp(1,1,ind_nu),n)
        enddo  !! over ind_nu
        deallocate(ps3,t,aa)
      enddo   !! over ind_k
      if(nproc_k/=1) call dgop(tmp,2*n*n*nd,'  +',comm_k)
c --------- Symmetrization --------------------------------------
      allocate(tmp1(n,n))
      do ii=1,num_sym_k(iq)
        ig=list_sym_k(ii,iq)
        do ind_nu=1,nd
          tmp1=tmp(:,:,ind_nu)
          call sym_chi_q(ig,tmp1,n,iq)
          chi0(:,:,ind_nu)=chi0(:,:,ind_nu)+tmp1
        enddo
      enddo
      deallocate(tmp1)
      chi0=chi0/num_sym_k(iq)
      deallocate(fifim,gb,tmp,zb02,ev02,zb2,ev2,indk0,wk0,ndk,nsk,zb1,
     &           ev1)
      if(isp1/=isp2) deallocate(zb01,ev01)
      call timel('**** CHI0_LDA_DYNAMIC_1Q finished **')
      end
