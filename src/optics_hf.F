      subroutine optics_hf(nrax,wrx,e)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nrax
      complex*16, intent(in) :: wrx(0:nrax)
      complex*16, intent(out) :: e(3,3,0:nrax)
	  integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,k0,inu,ibnd,npw,iq,
     &           iq0,npb,m,ne1,it,itr,m1,kq,kq0,nkq,neq
      real*8 :: de,fermi_dirac,fk,pi8,gtild(3,2),eq(3),sq2,v(3),df,ch,
     &          check_hermit
      complex*16 :: r011,r11,dr11,e11(6),ff,p11
      integer, allocatable :: ne(:)
      real*8, allocatable :: fifim(:,:,:,:,:),ve(:,:)
      complex*16, allocatable :: ea(:,:,:),vie(:,:,:),ppm(:,:,:),
     &                           ppe(:,:,:),aa(:,:,:,:),
     &                           zz(:,:,:,:),thet(:,:),ii(:,:,:,:),
     &                           jj(:,:,:,:),ck(:,:,:,:),ll(:,:,:),
     &                           r0m1(:),k0bb(:,:,:,:),tmp(:,:),
     &                           rm1(:),kbb(:,:,:,:),xx(:,:),drm1(:),
     &                           dkbb(:,:,:,:),fv(:),fxc(:),
     &                           thet_m(:,:,:),viet(:,:,:)
c ---------------- Z and A coefficients --------------------------      
      allocate(zz(nfun,nbndf_bnd,nqdiv,nspin))
      allocate(aa(nbasmpw,nbndf_bnd,nqdiv,nspin))
      zz=(0.d0,0.d0)
      aa=(0.d0,0.d0)
	  do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do ibnd=1,n_low_bnd(k,ispin)
            ib=ind_bands_bnd(ibnd,k,ispin)
	        zz(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
	        aa(:,ibnd,k,ispin)=ev_bnd(:,ib,ind_k,ispin)
          enddo
        enddo
      enddo
      if(nproc_k/=1) then
        call dgop(aa,2*nbasmpw*nbndf_bnd*npnt*nspin,'  +',comm_k)
        call dgop(zz,2*nfun*nbndf_bnd*npnt*nspin,'  +',comm_k)
      endif
      do ispin=1,nspin
        do k=npnt+1,nqdiv
	    k0=i_kref(k)
          call sym_z(zz(1,1,k,ispin),zz(1,1,k0,ispin),
     &               k_group(k),pnt(1,k),nfun,n_low_bnd(k0,ispin),
     &               maxb,io_lem,lmb,indbasa,limlb,0)
	      call sym_a(aa(1,1,k,ispin),aa(1,1,k0,ispin),k,k_group(k),
     &               pnt(1,k),nbasmpw,n_low_bnd(k0,ispin),nbask(k0),
     &               indgb(1,k),iplf_bk(1,k0))
        enddo
      enddo
c ----- Coulomb Interaction in OPB -------------------------------
      allocate(ne(nqdiv))
      allocate(ve(n_pbtot,nqdiv))
      allocate(vie(n_pbtot,n_pbtot,nqdiv))
      allocate(viet(n_pbtot,n_pbtot,nqdiv))
      do iq=1,nqdiv
        iq0=i_kref(iq)
        npw=nplwgw(iq0)
        n=n_pbmt+npw
c        if(iq==1) call v_coul_q0(n,npw,ve(1,iq),vie(1,1,iq),
c     &                           viet(1,1,iq),ne(1),nqdiv)
c        if(iq/=1) call v_coul_q(iq,ve(1,iq),vie(1,1,iq),viet(1,1,iq),
c     &                          ne(iq))
      enddo
      ne1=ne(1)-1
c ------- DFT Exchange-Correlation kernel ------------------------
      npw=nplwgw(1)
      n=n_pbmt+npw
      call resp_lda_prepare('100',.false.)
      allocate(thet_m(n,n,n_ixc_0))
      call theta_q_ss(1,n,thet_m)
      allocate(tmp(n,ne1))
      call zgemm('n','n',n,ne1,n,(1.d0,0.d0),thet_m,n,viet(1,2,1),
     &           n_pbtot,(0.d0,0.d0),tmp,n)
      deallocate(thet_m)
      allocate(thet(ne1,ne1))
      thet=(0.d0,0.d0)
      call zgemm('c','n',ne1,ne1,n,(1.d0,0.d0),viet(1,2,1),n_pbtot,tmp,
     &           n,(0.d0,0.d0),thet,ne1)
      deallocate(tmp)
c ----------------------------------------------------------------
      pi8=8.d0*pi
      allocate(fifim(maxel,maxel,n_pbmtm,nsort,nspin))
      do ispin=1,nspin
        call fi0_full(fifim(1,1,1,1,ispin),ispin,ispin)
      enddo
c ----------------------------------------------------------------
      allocate(ii(nbndf_bnd,nbndf_bnd,npnt,nspin))
      allocate(jj(nbndf_bnd,nbndf_bnd,npnt,nspin))
      allocate(ck(nbndf_bnd,nbndf_bnd,nqdiv,nspin))
      allocate(ll(nbndf_bnd,nqdiv,nspin))
      allocate(r0m1(ne1))
      allocate(k0bb(nbndf_bnd,nbndf_bnd,nqdiv,nspin))
      allocate(rm1(ne1))
      allocate(kbb(nbndf_bnd,nbndf_bnd,nqdiv,nspin))
      allocate(drm1(ne1))
      allocate(dkbb(nbndf_bnd,nbndf_bnd,nqdiv,nspin))
      allocate(xx(nbndf_bnd,nbndf_bnd))
      allocate(tmp(nbndf_bnd,nbndf_bnd))
      allocate(fv(ne1))
      allocate(fxc(ne1))
      sq2=sqrt(2.d0)/2.d0
      do inu=0,nrax
        do ispin=1,nspin
          do k=1,npnt
            do j=1,n_low_bnd(k,ispin)
              jb=ind_bands_bnd(j,k,ispin)
              de=e_bnd(jb,k,ispin)-chem_pot
              fk=fermi_dirac(de)
              do i=1,n_low_bnd(k,ispin)
                ib=ind_bands_bnd(i,k,ispin)
                de=e_bnd(ib,k,ispin)-chem_pot
                df=fk-fermi_dirac(de)
                de=e_bnd(jb,k,ispin)-e_bnd(ib,k,ispin)
                ii(i,j,k,ispin)=df/(wrx(inu)+de)
                jj(i,j,k,ispin)=df/(de*de-wrx(inu)**2)
              enddo
            enddo
          enddo
        enddo
c ---------------------------------------------------------------
        do it=1,6
          if(it==1) eq=(/1.d0,0.d0,0.d0/)
          if(it==2) eq=(/0.d0,1.d0,0.d0/)
          if(it==3) eq=(/0.d0,0.d0,1.d0/)
          if(it==4) eq=(/sq2,sq2,0.d0/)
          if(it==5) eq=(/sq2,0.d0,sq2/)
          if(it==6) eq=(/0.d0,sq2,sq2/)
          ck=(0.d0,0.d0)
          do ispin=1,nspin
            do k=1,nqdiv
              k0=i_kref(k)
              n=n_low_bnd(k0,ispin)
              allocate(ea(n,n,3))
              call bands_gradient(nbndf_bnd,n,k,ispin,zz(1,1,k,ispin),
     &                            aa(1,1,k,ispin),ea)
              ea=(0.d0,-2.d0)*ea
              do i=1,3
                ch=check_hermit(ea(1,1,i),n,n)
                ck(1:n,1:n,k,ispin)=ck(1:n,1:n,k,ispin)+eq(i)*ea(:,:,i)
              enddo
              deallocate(ea)
              do j=1,n_low_bnd(k0,ispin)
                jb=ind_bands_bnd(j,k0,ispin)
                de=e_bnd(jb,k0,ispin)-chem_pot
                fk=fermi_dirac(de)
                ll(j,k,ispin)=betta_t*fk*(1.d0-fk)*ck(j,j,k,ispin)
              enddo
            enddo
          enddo
c ----------- Initial values --------------------------------------
          r011=(0.d0,0.d0)
          r0m1=(0.d0,0.d0)
          k0bb=(0.d0,0.d0)
          do ispin=1,nspin
            do k=1,nqdiv
              k0=i_kref(k)
              n=n_low_bnd(k0,ispin)
              gtild=0.d0
              npb=n_pbmt+nplwgw(1)
              allocate(ppm(n,n,npb))
              call ppm_factors(k,k,1,fifim,zz(1,1,k,ispin),
     &                         zz(1,1,k,ispin),aa(1,1,k,ispin),
     &                         aa(1,1,k,ispin),ppm,gtild,npb,n,n,0,
     &                         nbndf_bnd)
              allocate(ppe(n,n,ne1))
	          call zgemm('n','n',n*n,ne1,npb,(1.d0,0.d0),ppm,n*n,
     &                   vie(1,2,1),n_pbtot,(0.d0,0.d0),ppe,n*n)
              deallocate(ppm)
              do j=1,n_low_bnd(k0,ispin)
                jb=ind_bands_bnd(j,k0,ispin)
                jbl=ind_block_k(jb,k0,ispin)
                do i=1,n_low_bnd(k0,ispin)
                  ib=ind_bands_bnd(i,k0,ispin)
                  ibl=ind_block_k(ib,k0,ispin)
                  de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                  if(ib==jb) then
                    k0bb(i,i,k,ispin)=ll(i,k,ispin)/wrx(inu)
                    r011=r011+ll(i,k,ispin)*ck(i,i,k,ispin)/wrx(inu)**2
                  else if(ibl/=jbl) then
                    k0bb(i,j,k,ispin)=-ck(i,j,k,ispin)*ii(i,j,k0,ispin)
     &                                /de
                    r011=r011+jj(i,j,k0,ispin)*abs(ck(i,j,k,ispin))**2
     &                                        /de
                    do m=1,ne1
                      r0m1(m)=r0m1(m)-conjg(ppe(i,j,m))*jj(i,j,k0,ispin)
     &                               *ck(i,j,k,ispin)
                    enddo
                  endif
                enddo
              enddo
              deallocate(ppe)
            enddo   !! over k
          enddo   !! over ispin
          k0bb=k0bb/sqrt(amega)
          r011=r011/nqdiv/amega
          r0m1=r0m1/nqdiv/sqrt(amega)
          if(nspin==1) then
            r011=2.d0*r011
            r0m1=2.d0*r0m1
          endif
c          p11=r011
c -----------------------------------------------------------------
          r11=r011
          rm1=r0m1
          kbb=k0bb
          do itr=1,iter_ladder_xi
            dr11=(0.d0,0.d0)
            drm1=(0.d0,0.d0)
            dkbb=(0.d0,0.d0)
            fv=(0.d0,0.d0)
            fxc=(0.d0,0.d0)
            do m=1,ne1
              do m1=1,ne1
                fxc(m)=fxc(m)+thet(m,m1)*rm1(m1)
              enddo
              fv(m)=ve(1+m,1)*rm1(m)
            enddo
            do ispin=1,nspin
              do k=1,nqdiv
                k0=i_kref(k)
                n=n_low_bnd(k0,ispin)
                xx=(0.d0,0.d0)
                do iq=1,nqdiv
                  iq0=i_kref(iq)
                  v=pnt(:,k)-pnt(:,iq)
                  call zone1_number(v,rb0,ndiv,kq)
                  kq=index_k1(kq)    !! for K-Q
                  kq0=i_kref(kq)
                  nkq=n_low_bnd(kq0,ispin)
                  gtild(:,1)=0.d0
                  gtild(:,2)=v-pnt(:,kq)
                  npb=n_pbmt+nplwgw(iq0)
                  allocate(ppm(n,nkq,npb))
                  call ppm_factors(k,kq,iq,fifim,zz(1,1,k,ispin),
     &                             zz(1,1,kq,ispin),aa(1,1,k,ispin),
     &                             aa(1,1,kq,ispin),ppm,gtild,npb,n,nkq,
     &                             0,nbndf_bnd)
                  neq=ne(iq)
                  allocate(ppe(n,nkq,neq))
	              call zgemm('n','n',n*nkq,neq,npb,(1.d0,0.d0),ppm,
     &                       n*nkq,vie(1,1,iq),n_pbtot,(0.d0,0.d0),ppe,
     &                       n*nkq)
                  deallocate(ppm)
                  do m=1,neq
	              call zgemm('n','c',nkq,n,nkq,(1.d0,0.d0),
     &                         kbb(1,1,kq,ispin),nbndf_bnd,ppe(1,1,m),
     &                         n,(0.d0,0.d0),tmp,nbndf_bnd)
	              call zgemm('n','n',n,n,nkq,ve(m,iq),ppe(1,1,m),n,
     &                         tmp,nbndf_bnd,(1.d0,0.d0),xx,nbndf_bnd)
                  enddo
                  deallocate(ppe)
                enddo   !! over iq
                xx=alpha_hbr*xx/nqdiv
c -------------------------------------------------------------------
                npb=n_pbmt+nplwgw(1)
                allocate(ppm(n,n,npb))
                gtild=0.d0
                call ppm_factors(k,k,1,fifim,zz(1,1,k,ispin),
     &                           zz(1,1,k,ispin),aa(1,1,k,ispin),
     &                           aa(1,1,k,ispin),ppm,gtild,npb,n,n,0,
     &                           nbndf_bnd)
                allocate(ppe(n,n,ne1))
	            call zgemm('n','n',n*n,ne1,npb,(1.d0,0.d0),ppm,n*n,
     &                     vie(1,2,1),n_pbtot,(0.d0,0.d0),ppe,n*n)
                deallocate(ppm)
c -------------------------------------------------------------------
                do j=1,n_low_bnd(k0,ispin)
                  jb=ind_bands_bnd(j,k0,ispin)
                  jbl=ind_block_k(jb,k0,ispin)
                  do i=1,n_low_bnd(k0,ispin)
                    ib=ind_bands_bnd(i,k0,ispin)
                    ibl=ind_block_k(ib,k0,ispin)
                    de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                    if(ib==jb) then
                      dr11=dr11-ll(i,k,ispin)*xx(i,i)/wrx(inu)
     &                         /sqrt(amega)/nqdiv
                    else if(ibl/=jbl) then
                      ff=(0.d0,0.d0)
                      do m=1,ne1
                        ff=ff+ppe(i,j,m)*fxc(m)
                      enddo
                      dkbb(i,j,k,ispin)=ii(i,j,k0,ispin)*(ff+xx(i,j))
                      dr11=dr11+jj(i,j,k0,ispin)*conjg(ck(i,j,k,ispin))
     &                         *ff/sqrt(amega)/nqdiv
                      do m=1,ne1
                        drm1(m)=drm1(m)-conjg(ppe(i,j,m))
     &                         *jj(i,j,k0,ispin)*de*ff/nqdiv
                      enddo
                      dr11=dr11+ii(i,j,k0,ispin)*xx(i,j)/de/sqrt(amega)
     &                         /nqdiv*conjg(ck(i,j,k,ispin))
                      do m=1,ne1
                        drm1(m)=drm1(m)-conjg(ppe(i,j,m))
     &                         *ii(i,j,k0,ispin)*xx(i,j)/nqdiv
                      enddo
                    endif
                  enddo   !! over i
                enddo   !! over j
                deallocate(ppe)
              enddo  !! over k
            enddo  !! over ispin
            if(nspin==1) then
              dr11=2.d0*dr11
              drm1=2.d0*drm1
            endif
            kbb=k0bb+dkbb
            r11=r011+dr11
            rm1=r0m1+drm1
            print *,inu,itr,r11
          enddo   !! over itr
          p11=r11   !p11+r11      
          e11(it)=(1.d0,0.d0)-8.d0*pi*p11
          e11(it)=(1.d0,0.d0)/e11(it)
        enddo  !! over it
        e(1,1,inu)=(1.d0,0.d0)/e11(1)
        e(2,2,inu)=(1.d0,0.d0)/e11(2)
        e(3,3,inu)=(1.d0,0.d0)/e11(3)
        e(1,2,inu)=(1.d0,0.d0)/e11(4)-0.5d0*(e(1,1,inu)+e(2,2,inu))
        e(1,3,inu)=(1.d0,0.d0)/e11(5)-0.5d0*(e(1,1,inu)+e(3,3,inu))
        e(2,3,inu)=(1.d0,0.d0)/e11(6)-0.5d0*(e(2,2,inu)+e(3,3,inu))
        e(2,1,inu)=e(1,2,inu)
        e(3,1,inu)=e(1,3,inu)
        e(3,2,inu)=e(2,3,inu)
      enddo  !! over inu
      deallocate(aa,zz,ne,ve,vie,thet,fifim,ii,jj,ck,ll,r0m1,k0bb,rm1,
     &           kbb,xx,tmp,drm1,dkbb,fv,fxc,viet)
      end
