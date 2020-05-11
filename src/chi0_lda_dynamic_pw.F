      subroutine chi0_lda_dynamic_pw(ch,wrx,xi0)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      real*8, intent(in) :: wrx(0:nrax_chi)
      complex*16, intent(inout) :: xi0(0:nrax_chi)
      integer :: k,ind_k,i,j,n1,n2,kq,kq0,k0,ks,i_nu,nq0,isp1,isp2,iq,
     &           ii
      real*8 :: gtild(3),de,v(3),fermi_dirac,d,w,ek,ekq,fk,fkq,d0
      complex*16 :: nuc
      integer, allocatable :: indk0(:),ndk(:),nsk(:)
      real*8, allocatable :: gb(:,:),fifi_j(:,:,:,:),wk0(:)
      complex*16, allocatable :: ps3(:,:),zb02(:,:,:),ev02(:,:,:),
     &                           zb2(:,:),ev2(:,:),zb01(:,:,:),
     &                           ev01(:,:,:),zb1(:,:),ev1(:,:),
     &                           ff(:,:,:),gre(:,:,:),ea(:,:,:)
      isp1=1
      isp2=1
      if(ch=='___XY') isp1=2
      iq=q_suscept
      allocate(fifi_j(maxlfun,maxlfun,0:2*maxb,nsort))
      call fifi_j_prepare_0(fifi_j,isp1,isp2,pnt(1,iq))
      allocate(ff(maxel,maxel,natom))
      call psi_psi_pw_full(fifi_j,ff,pnt(1,iq))
      deallocate(fifi_j)
      allocate(gb(nbndf,2))
      xi0=(0.d0,0.d0)
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
      do ind_k=1,ndk(me_k+1)
        ks=nsk(me_k+1)+ind_k
        k=indk0(ks)
        k0=i_kref(k)
        v=pnt(:,k)-pnt(:,iq)
        call zone1_number(v,rb0,ndiv,kq)
        kq=index_k1(kq)    !! for K-Q
        kq0=i_kref(kq)
        gtild=v-pnt(:,kq)
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
        allocate(ps3(n1,n2))
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
        call psi_psi_pw_end(n1,n2,gtild,ev1,zb1,ev2,zb2,k,kq,ff,ps3,1)
        if(nn_k_int==1) then
          do j=1,n2
            do i=1,n1
              de=e_bnd(j,kq0,isp2)-e_bnd(i,k0,isp1)
              d=real(ps3(i,j)*conjg(ps3(i,j)))*(gb(j,2)-gb(i,1))*wk0(ks)
              do i_nu=0,nrax_chi
                nuc=dcmplx(wrx(i_nu),e_small_bos)
                xi0(i_nu)=xi0(i_nu)+d/(nuc+de)
              enddo
            enddo
          enddo
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
          do j=1,n2
            do i=1,n1
              d0=real(ps3(i,j)*conjg(ps3(i,j)))*w
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
                d=d0*(fkq-fk)
                do i_nu=0,nrax_chi
                  nuc=dcmplx(wrx(i_nu),e_small_bos)
                  xi0(i_nu)=xi0(i_nu)+d/(nuc+de)
                enddo
              enddo
            enddo
          enddo
          deallocate(gre)
        endif
        deallocate(ps3)
      enddo   !! over ind_k
      if(nproc_k/=1) call dgop(xi0,2*(nrax_chi+1),'  +',comm_k)
      deallocate(ff,gb,zb02,ev02,zb2,ev2,indk0,wk0,ndk,nsk,zb1,ev1)
      if(isp1/=isp2) deallocate(zb01,ev01)
      if(ch=='___00') xi0=2.d0*xi0
      if(ch=='___ZZ') xi0=2.d0*xi0
      if(ch=='___XY') xi0=4.d0*xi0
      call timel('**** CHI0_LDA_DYNAMIC_PW finished **')
      end
