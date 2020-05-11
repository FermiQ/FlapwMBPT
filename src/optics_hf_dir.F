      subroutine optics_hf_dir(nrax,wrx,e,mode)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      character*1, intent(in) :: mode
      integer, intent(in) :: nrax
      complex*16, intent(in) :: wrx(0:nrax)
      complex*16, intent(out) :: e(3,3,0:nrax)
	  integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,k0,inu,npw,iq,k0a,
     &           iq0,npb,m,ne1,it,neq,nmm,info,jspin,ind_nu,ka,k10a,
     &           k1,k10,n1,j1,i1,jnd,ind,k1a,iq0a,iqa,kc
      real*8 :: de,fermi_dirac,fk,pi8,gtild(3,2),eq(3),sq2,v(3),df,ch,
     &          check_hermit
      complex*16 :: e11(6),ff,p11,dp11
      integer, allocatable :: ne(:),indm(:,:,:,:),ipiv(:),nd(:),ns(:)
      real*8, allocatable :: fifim(:,:,:,:,:),ve(:,:)
      complex*16, allocatable :: ea(:,:,:),vie(:,:,:),ppm(:,:,:),
     &                           ppe(:,:,:),aa(:,:,:,:),
     &                           zz(:,:,:,:),thet(:,:,:,:),ii(:,:,:,:),
     &                           jj(:,:,:,:),ck(:,:,:,:),ll(:,:,:),
     &                           thet_m(:,:,:,:),tmp(:,:),
     &                           viet(:,:,:),psi(:,:,:,:),
     &                           aaa(:,:),bb(:),vv0(:,:,:,:),
     &                           fxc(:,:),pp(:,:,:,:,:),ttt(:,:,:),
     &                           pab(:,:,:)
      allocate(indm(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      it=0
      do ispin=1,nspin
        do k=1,nqdiv_c
          k0=i_kref_c(k)
          k0a=k_a_from_c(k0)
          n=n_low_bnd(k0a,ispin)
          do j=1,n
            do i=1,n
              it=it+1
              indm(i,j,k,ispin)=it
            enddo
          enddo
        enddo
      enddo
      nmm=it
      allocate(nd(nproc))
      allocate(ns(nproc))
      call size_shift_par(nrax+1,nproc,nd,ns)
c ---------------- Z and A coefficients --------------------------      
      allocate(zz(nfun,nbndf_bnd,nqdiv_c,nspin))
      allocate(aa(nbasmpw,nbndf_bnd,nqdiv_c,nspin))  
      zz=(0.d0,0.d0)
      aa=(0.d0,0.d0)
	  do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          kc=k_c_from_a(k)
          if(kc==0) cycle
          do i=1,n_low_bnd(k,ispin)
            ib=ind_bands_bnd(i,k,ispin)
	        zz(:,i,kc,ispin)=z_bnd(:,ib,ind_k,ispin)
	        aa(:,i,kc,ispin)=ev_bnd(:,ib,ind_k,ispin)
          enddo
        enddo
        if(nproc_k/=1) then
          call dgop(aa(1,1,1,ispin),2*nbasmpw*nbndf_bnd*npnt_c,'  +',
     &              comm_k)
          call dgop(zz(1,1,1,ispin),2*nfun*nbndf_bnd*npnt_c,'  +',
     &              comm_k)
        endif
        do k=npnt_c+1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k)
          k0a=k_a_from_c(k0)
          call sym_z(zz(1,1,k,ispin),zz(1,1,k0,ispin),k_group_c(k),
     &               pnt_c(1,k),nfun,n_low_bnd(k0a,ispin),maxb,io_lem,
     &               lmb,indbasa,limlb,0)
	      call sym_a(aa(1,1,k,ispin),aa(1,1,k0,ispin),ka,k_group_c(k),
     &               pnt_c(1,k),nbasmpw,n_low_bnd(k0a,ispin),nbask(k0a),
     &               indgb(1,ka),iplf_bk(1,k0a))
        enddo
      enddo
      if(mode/='1') then
c ----- Coulomb Interaction in OPB -------------------------------
        allocate(ne(nqdiv_c))
        allocate(ve(n_pbtot,nqdiv_c))
        allocate(vie(n_pbtot,n_pbtot,nqdiv_c))
        allocate(viet(n_pbtot,n_pbtot,nqdiv_c))
        do iq=1,nqdiv_c
          iq0=i_kref_c(iq)
          k0a=k_a_from_c(iq0)
          ka=k_a_from_c(iq)
          npw=nplwgw(k0a)
          n=n_pbmt+npw
c          if(iq==1) call v_coul_q0(n,npw,ve(1,iq),vie(1,1,iq),
c     &                             viet(1,1,iq),ne(1),nqdiv_c)
c          if(iq/=1) call v_coul_q(ka,ve(1,iq),vie(1,1,iq),viet(1,1,iq),
c     &                            ne(iq))
        enddo
        ne1=ne(1)-1
c ------- DFT Exchange-Correlation kernel ------------------------
        npw=nplwgw(1)
        n=n_pbmt+npw
        n_ixc_0=4
        call theta_xc(1)
        allocate(thet_m(n,n,2,2))
        call theta_q_ss(1,n,thet_m)
        allocate(tmp(n,ne(1)))
        if(nspin==1) thet_m(:,:,1,1)=thet_m(:,:,1,1)+thet_m(:,:,1,2)
        allocate(thet(ne(1),ne(1),nspin,nspin))
        thet=(0.d0,0.d0)
        do jspin=1,nspin
          do ispin=1,nspin
            call zgemm('n','n',n,ne(1),n,(1.d0,0.d0),
     &                 thet_m(1,1,ispin,jspin),n,viet,n_pbtot,
     &                 (0.d0,0.d0),tmp,n)
            call zgemm('c','n',ne(1),ne(1),n,(1.d0,0.d0),viet,n_pbtot,
     &                 tmp,n,(0.d0,0.d0),thet(1,1,ispin,jspin),ne(1))
          enddo
        enddo
        deallocate(thet_m,tmp)
c ----------------------------------------------------------------
        pi8=8.d0*pi
        allocate(fifim(maxel,maxel,n_pbmtm,nsort,nspin))
        do ispin=1,nspin
          call fi0_full(fifim(1,1,1,1,ispin),ispin,ispin)
        enddo
      endif
c ----------------------------------------------------------------
      allocate(ii(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
      allocate(jj(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
      allocate(ck(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      allocate(ll(nbndf_bnd,nqdiv_c,nspin))
      allocate(psi(nbndf_bnd,nbndf_bnd,nqdiv_c,nspin))
      sq2=sqrt(2.d0)/2.d0
      if(mode/='1') then
c ----------- Matrix elements (4-point) of f_xc ---------------------
        allocate(pp(nbndf_bnd,nbndf_bnd,ne(1),nqdiv_c,nspin))
        do ispin=1,nspin
          do k=1,nqdiv_c
            k0=i_kref_c(k)
            k0a=k_a_from_c(k0)
            ka=k_a_from_c(k)
            n=n_low_bnd(k0a,ispin)
            gtild=0.d0
            npb=n_pbmt+nplwgw(1)
            allocate(ppm(n,n,npb))
            call ppm_factors(ka,ka,1,fifim,zz(1,1,k,ispin),
     &                       zz(1,1,k,ispin),aa(1,1,k,ispin),
     &                       aa(1,1,k,ispin),ppm,gtild,npb,n,n,0,
     &                       nbndf_bnd)
            allocate(ppe(n,n,ne(1)))
	        call zgemm('n','n',n*n,ne(1),npb,(1.d0,0.d0),ppm,n*n,
     &                 vie,n_pbtot,(0.d0,0.d0),ppe,n*n)
            deallocate(ppm)
            pp(1:n,1:n,:,k,ispin)=ppe
            deallocate(ppe)
          enddo
        enddo
        allocate(fxc(nmm,nmm))
        do jspin=1,nspin
          do ispin=1,nspin
            do k1=1,nqdiv_c
              k10=i_kref_c(k1)
              k10a=k_a_from_c(k10)
              n1=n_low_bnd(k10a,ispin)
              allocate(ppm(n1,n1,ne(1)))
              ppm=pp(1:n1,1:n1,:,k1,jspin)
              allocate(ttt(ne(1),n1,n1))
	          call zgemm('n','c',ne(1),n1*n1,ne(1),(1.d0,0.d0),
     &                   thet(1,1,ispin,jspin),ne(1),ppm,n1*n1,
     &                   (0.d0,0.d0),ttt,ne(1))
              do k=1,nqdiv_c
                k0=i_kref_c(k)
                k0a=k_a_from_c(k0)
                n=n_low_bnd(k0a,ispin)
                allocate(ppe(n,n,ne(1)))
                ppe=pp(1:n,1:n,:,k,ispin)
                allocate(vv0(n,n,n1,n1))
                call zgemm('n','n',n*n,n1*n1,ne(1),(1.d0,0.d0),ppe,n*n,
     &                     ttt,ne(1),(0.d0,0.d0),vv0,n*n)
                deallocate(ppe)
                do j1=1,n1
                  do i1=1,n1
                    jnd=indm(i1,j1,k1,jspin)
                    do j=1,n
                      do i=1,n
                        ind=indm(i,j,k,ispin)
                        fxc(ind,jnd)=vv0(i,j,i1,j1)
                      enddo
                    enddo
                  enddo
                enddo
                deallocate(vv0)
              enddo   !! over k
              deallocate(ttt,ppm)
            enddo   !! over k1
          enddo   !! over ispin
        enddo  !! over jspin
        deallocate(pp)
        if(mode=='3') then
c ----------- Matrix elements (4-point) of V --------------------- 
          do jspin=1,nspin
            do ispin=1,nspin
              if(ispin/=jspin) cycle
              do k1=1,nqdiv_c
                k10=i_kref_c(k1)
                k10a=k_a_from_c(k10)
                k1a=k_a_from_c(k1)
                n1=n_low_bnd(k10a,ispin)
                do k=1,nqdiv_c
                  k0=i_kref_c(k)
                  k0a=k_a_from_c(k0)
                  ka=k_a_from_c(k)
                  v=pnt_c(:,k)-pnt_c(:,k1)
                  call zone1_number(v,rb0_c,ndiv_c,iq)
                  iq=index_k1_c(iq)    !! for K-K'
                  iq0=i_kref_c(iq)
                  iq0a=k_a_from_c(iq0)
                  iqa=k_a_from_c(iq)
                  n=n_low_bnd(k0a,ispin)
                  gtild(:,1)=0.d0
                  gtild(:,2)=v-pnt_c(:,iq)
                  npb=n_pbmt+nplwgw(iq0a)
                  neq=ne(iq)
                  allocate(ppe(n,n1,neq))
                  allocate(ppm(n,n1,npb))
                  call ppm_factors(ka,k1a,iqa,fifim,zz(1,1,k,ispin),
     &                             zz(1,1,k1,ispin),aa(1,1,k,ispin),
     &                             aa(1,1,k1,ispin),ppm,gtild,npb,n,n1,
     &                             0,nbndf_bnd)
	              call zgemm('n','n',n*n1,neq,npb,(1.d0,0.d0),ppm,
     &                       n*n1,vie(1,1,iq),n_pbtot,(0.d0,0.d0),ppe,
     &                       n*n1)
                  deallocate(ppm)
                  if(iq==1) then
                    ppe(:,:,1)=(0.d0,0.d0)
                    do i=1,n
                      ppe(i,i,1)=(1.d0,0.d0)/sqrt(amega)
                    enddo
                  endif
                  allocate(ppm(n,n1,neq))
                  do m=1,neq
                    ppm(:,:,m)=ppe(:,:,m)*ve(m,iq)
                  enddo
                  allocate(vv0(n,n1,n,n1))
              call zgemm('n','c',n*n1,n*n1,neq,(1.d0,0.d0),ppe,n*n1,
     &                       ppm,n*n1,(0.d0,0.d0),vv0,n*n1)
                  deallocate(ppe,ppm)
                  do j1=1,n1
                    do i1=1,n1
                      jnd=indm(i1,j1,k1,jspin)
                      do j=1,n
                        do i=1,n
                          ind=indm(i,j,k,ispin)
                          fxc(ind,jnd)=fxc(ind,jnd)
     &                                -vv0(i,i1,j,j1)*alpha_hbr
                        enddo
                      enddo
                    enddo
                  enddo
                  deallocate(vv0)
                enddo   !! over k
              enddo   !! over k1
            enddo   !! over ispin
          enddo  !! over jspin
        endif
        fxc=fxc/nqdiv
        call timel('* Fxc+VV in OPTICS_HF_DIR finished *')
      endif
c ------- Head of non-interacting polarizability ----------------
      allocate(pab(3,3,0:nrax))
      call p0_head_tensor(nrax,wrx,pab)
      call timel('* Head of P0 finished **************')
c ---------------------------------------------------------------
      e=(0.d0,0.d0)
	  do ind_nu=1,nd(me+1)
	    inu=ns(me+1)+ind_nu-1
        do ispin=1,nspin
          do k=1,npnt_c
            ka=k_a_from_c(k)
            do j=1,n_low_bnd(ka,ispin)
              jb=ind_bands_bnd(j,ka,ispin)
              de=e_bnd(jb,ka,ispin)-chem_pot
              fk=fermi_dirac(de)
              do i=1,n_low_bnd(ka,ispin)
                ib=ind_bands_bnd(i,ka,ispin)
                de=e_bnd(ib,ka,ispin)-chem_pot
                df=fk-fermi_dirac(de)
                de=e_bnd(jb,ka,ispin)-e_bnd(ib,ka,ispin)
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
          psi=(0.d0,0.d0)
          do ispin=1,nspin
            do k=1,nqdiv_c
              k0=i_kref_c(k)
              k0a=k_a_from_c(k0)
              ka=k_a_from_c(k)
              n=n_low_bnd(k0a,ispin)
              allocate(ea(n,n,3))
              call bands_gradient(nbndf_bnd,n,ka,ispin,zz(1,1,k,ispin),
     &                            aa(1,1,k,ispin),ea)
              ea=(0.d0,-2.d0)*ea
              do i=1,3
                ch=check_hermit(ea(1,1,i),n,n)
                ck(1:n,1:n,k,ispin)=ck(1:n,1:n,k,ispin)+eq(i)*ea(:,:,i)
              enddo
              deallocate(ea)
              do j=1,n_low_bnd(k0a,ispin)
                jb=ind_bands_bnd(j,k0a,ispin)
                de=e_bnd(jb,k0a,ispin)-chem_pot
                fk=fermi_dirac(de)
                ll(j,k,ispin)=betta_t*fk*(1.d0-fk)*ck(j,j,k,ispin)
              enddo
              do j=1,n_low_bnd(k0a,ispin)
                jb=ind_bands_bnd(j,k0a,ispin)
                jbl=ind_block_k(jb,k0a,ispin)
                do i=1,n_low_bnd(k0a,ispin)
                  ib=ind_bands_bnd(i,k0a,ispin)
                  ibl=ind_block_k(ib,k0a,ispin)
                  de=e_bnd(jb,k0a,ispin)-e_bnd(ib,k0a,ispin)
                  if(ib==jb) then
                    if(j==i) psi(i,j,k,ispin)=ll(i,k,ispin)/wrx(inu)
     &                                                     /sqrt(amega)
                  else if(ibl/=jbl) then
                    psi(i,j,k,ispin)=-ck(i,j,k,ispin)*ii(i,j,k0,ispin)
     &                               /de/sqrt(amega)
                  endif
                enddo
              enddo
            enddo
          enddo
c ----------- Initial values --------------------------------------
          p11=(0.d0,0.d0)
          do j=1,3
            do i=1,3
              p11=p11+eq(i)*eq(j)*pab(i,j,inu)
            enddo
          enddo
c -----------------------------------------------------------------
          if(mode/='1') then
            allocate(aaa(nmm,nmm))
            allocate(bb(nmm))
            allocate(ipiv(nmm))
            bb=(0.d0,0.d0)
            aaa=(0.d0,0.d0)
            do jspin=1,nspin
              do ispin=1,nspin
                do k1=1,nqdiv_c
                  k10=i_kref_c(k1)
                  k10a=k_a_from_c(k10)
                  n1=n_low_bnd(k10a,ispin)
                  do k=1,nqdiv_c
                    k0=i_kref_c(k)
                    k0a=k_a_from_c(k0)
                    n=n_low_bnd(k0a,ispin)
                    do j1=1,n1
                      do i1=1,n1
                        jnd=indm(i1,j1,k1,jspin)
                        do j=1,n
                          do i=1,n
                            ind=indm(i,j,k,ispin)
                            if(ind==jnd) aaa(ind,jnd)=aaa(ind,jnd)
     &                                               +(1.d0,0.d0)
                            ff=fxc(ind,jnd)
                            aaa(ind,jnd)=aaa(ind,jnd)
     &                                  +ff*ii(i1,j1,k10,jspin)
                            bb(ind)=bb(ind)-ff*psi(i1,j1,k1,jspin)
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo  !! over k1
              enddo  !! over ispin
            enddo   !! over jspin
            call zgesv(nmm,1,aaa,nmm,ipiv,bb,nmm,info)
            deallocate(aaa,ipiv)
            dp11=(0.d0,0.d0)
            do ispin=1,nspin
              do k=1,nqdiv_c
                k0=i_kref_c(k)
                k0a=k_a_from_c(k0)
                n=n_low_bnd(k0a,ispin)
                do j=1,n
                  jb=ind_bands_bnd(j,k0a,ispin)
                  jbl=ind_block_k(jb,k0a,ispin)
                  do i=1,n
                    ib=ind_bands_bnd(i,k0a,ispin)
                    ibl=ind_block_k(ib,k0a,ispin)
                    de=e_bnd(jb,k0a,ispin)-e_bnd(ib,k0a,ispin)
                    ind=indm(i,j,k,ispin)
                    ff=(0.d0,0.d0)
                    if(i==j) then
                      ff=ll(i,k,ispin)/wrx(inu)
                    else if(ibl/=jbl) then
                      ff=-conjg(ck(i,j,k,ispin))*ii(i,j,k0,ispin)/de
                    endif
                    dp11=dp11-ff*bb(ind)/nqdiv_c/sqrt(amega)
                  enddo
                enddo
              enddo   !! over k
            enddo   !! over ispin
            deallocate(bb)
            if(nspin==1) dp11=2.d0*dp11
            p11=p11+dp11
          endif
          e11(it)=(1.d0,0.d0)-8.d0*pi*p11
          e11(it)=(1.d0,0.d0)/e11(it)
c          print *,inu,p11
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
      enddo  !! over ind_nu
#ifdef MPI
      if(nproc/=1) call dgop(e,18*(nrax+1),'  +',MPI_COMM_WORLD)
#endif      
      deallocate(aa,zz,ii,jj,ck,ll,indm,psi,nd,ns,pab)
      if(mode/='1') deallocate(thet,ne,ve,vie,viet,fifim,fxc)
      call timel('***** OPTICS_HF_DIR finished *******')
      end
