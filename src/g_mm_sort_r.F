      subroutine g_mm_sort_r(gr,it,ind_tau,ispin,nf1,nf2,nd_nrr,
     &                       n0_nrr,nrr_reds,nat1,nat2,isort,jsort,
     &                       nfun1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: it,ispin,ind_tau,nf1,nf2,nd_nrr,n0_nrr,
     &                       nat1,nat2,nrr_reds(3,nqdiv*nat1*nat2),
     &                       isort,jsort,nfun1
      real*8, intent(out) :: gr(nf1,nf2,nd_nrr)
      integer :: k0,nbnd,ind,j,ib,ind_k,ig,i_tau,j0,ir,ir0,ieapw,
     &           i0,i1,ind_ir0,i,m,nls,l,npr,ndim0,dli,dli1,
     &           nd1,nd,i3,i2,ind_ls,ls,l1,li1,ne1,ns1,indls1,li,ie,
     &           ie00,ie11,iend,jend,ii,jatom1,iatom1,
     &           ie1,ii1,ir1,jatom0,jat0,iatom0,iat0,m1,lm1,
     &           km1,lm,km,nz1,iat1,jat1
      real*8 :: pi2,tu,tt(3),phase,v0(3),v(3),v1(3)
      integer, allocatable :: adr_ls(:),nd_nls(:),nm_nls(:,:)
      real*8, allocatable :: gx_tau(:),g_tmp(:,:,:,:),trf(:,:,:),
     &                       g_mm(:,:,:),aa(:,:),bb(:,:)
      complex*16, allocatable :: tmp1(:,:),tm(:,:),tmp2(:,:),
     &                           z_ord(:,:,:)
      allocate(tmp1(nbndf,nbndf))
      pi2=pi+pi
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      nz1=io_lem(iat_1(isort))
c -----------------------------------------------------------------
      allocate(adr_ls(lmb(jsort)+1))
      ind=0
      do l=0,lmb(jsort)
        ind=ind+1
        adr_ls(ind)=l
      enddo
      nls=ind
      npr=nproc_b
      ndim0=nls/npr
      if(ndim0*npr<nls) ndim0=ndim0+1
      allocate(nd_nls(npr))
      allocate(nm_nls(ndim0,npr))
      call size_shift_par_load(nls,npr,nd_nls,nm_nls,ndim0)
      allocate(trf(ndim3_k(me_k+1),nqdiv,2))
      ir=0
      do i3=0,ndiv(3)-1
        do i2=0,ndiv(2)-1
          do i1=0,ndiv(1)-1
            ir=ir+1
            tt=rbas(:,1)*i1+rbas(:,2)*i2+rbas(:,3)*i3
            do ind_k=1,ndim3_k(me_k+1)
              k0=n3_mpi_k(me_k+1)+ind_k
              phase=pi2*dot_product(pnt(:,k0),tt)
              trf(ind_k,ir,1)=wgt(k0)*cos(phase)
              trf(ind_k,ir,2)=wgt(k0)*sin(phase)
            enddo
          enddo
        enddo
      enddo
c -----------------------------------------------------------------
      gr=0.d0
      allocate(gx_tau(nbndf))
      allocate(z_ord(nfun,nbndf,ndim3_k(me_k+1)))
      do ind_k=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_k
        nbnd=n_bnd(k0,ispin)
        do j=1,nfun
          j0=nfun_order(j)
          do ib=1,nbnd
            z_ord(j,ib,ind_k)=z_bnd(j0,ib,ind_k,ispin)
          enddo
        enddo
      enddo
c -----------------------------------------------------------------
      do ind_ls=1,nd_nls(me_b+1)
        ls=nm_nls(ind_ls,me_b+1)
        l1=adr_ls(ls)
        li1=l1+1
        nd1=l1+l1+1
        ne1=ntle(l1,jsort)
        ieapw=n_e_apw(l1,jsort)
        if(augm(ieapw,l1,jsort)=='LAW'.or.
     &     augm(ieapw,l1,jsort)=='AWL') ne1=ne1+1
        ns1=nd1*ne1*nat2
        indls1=nfun_ls_adr(li1,jsort)
        allocate(g_tmp(nfun1,ns1,ndim3_k(me_k+1),2))
        allocate(tmp2(nbndf,ns1))
        allocate(tm(nfun1,ns1))
        do ind_k=1,ndim3_k(me_k+1)
          k0=n3_mpi_k(me_k+1)+ind_k
          nbnd=n_bnd(k0,ispin)
          if(ubi=='dft'.or.ubi==' hf') then
            tu=tau_mesh(i_tau)
            call g_x_tau(ispin,k0,gx_tau,tu,nbnd,chem_pot)
            do j=1,ns1
              j0=indls1+j-1
              do ib=1,nbnd
                tmp2(ib,j)=gx_tau(ib)*conjg(z_ord(j0,ib,ind_k))
              enddo
            enddo
          else
            call ferm_unpack_tau(tmp1,g_full(1,1,1,ind_tau,ind_k,ispin),
     &                           nbnd,nbndf,nbndf,it)
            call zgemm('n','c',nbnd,ns1,nbnd,(1.d0,0.d0),tmp1,nbndf,
     &                 z_ord(indls1,1,ind_k),nfun,(0.d0,0.d0),tmp2,
     &                 nbndf)
          endif
          call zgemm('n','n',nfun1,ns1,nbnd,(1.d0,0.d0),
     &               z_ord(nz1,1,ind_k),nfun,tmp2,nbndf,(0.d0,0.d0),tm,
     &               nfun1)
          g_tmp(:,:,ind_k,1)=real(tm)
          g_tmp(:,:,ind_k,2)=imag(tm)
        enddo   !! over ind_k
        deallocate(tmp2,tm)
        allocate(g_mm(nfun1,ns1,nqdiv))
        call dgemm('n','n',nfun1*ns1,nqdiv,ndim3_k(me_k+1),1.d0,
     &             g_tmp(1,1,1,1),nfun1*ns1,trf(1,1,1),ndim3_k(me_k+1),
     &             0.d0,g_mm,nfun1*ns1)
        call dgemm('n','n',nfun1*ns1,nqdiv,ndim3_k(me_k+1),-1.d0,
     &             g_tmp(1,1,1,2),nfun1*ns1,trf(1,1,2),ndim3_k(me_k+1),
     &             1.d0,g_mm,nfun1*ns1)
        deallocate(g_tmp)
        if(nproc_k/=1) call dgop(g_mm,nfun1*ns1*nqdiv,'  +',comm_k)
c ------------------------------------------------------------------
        dli1=l1*(nd1-2)*nd1/3+1
        do ind_ir0=1,nd_nrr
          ir0=n0_nrr+ind_ir0
          jat1=nrr_reds(3,ir0)
          jatom1=iat_sort(jat1,jsort)
          ir=nrr_reds(1,ir0)
          v0=rr_3(:,ir)
          iat1=nrr_reds(2,ir0)
          iatom1=iat_sort(iat1,isort)
          ie1=0
          do ie11=1,ntle(l1,jsort)
            jend=1
            if(augm(ie11,l1,jsort)=='LAW') jend=2
            if(augm(ie11,l1,jsort)=='AWL') jend=2
            do ii1=1,jend
              ie1=ie1+1
              do l=0,lmb(isort)
                li=l+1
                nd=l+l+1
                allocate(aa(nd,nd1))
                allocate(bb(nd,nd1))
                dli=l*(nd-2)*nd/3+1
                ie=0
                do ie00=1,ntle(l,isort)
                  iend=1
                  if(augm(ie00,l,isort)=='LAW') iend=2
                  if(augm(ie00,l,isort)=='AWL') iend=2
                  do ii=1,iend
                    ie=ie+1
                    bb=0.d0
                    do ig=1,ngroup
                      v=v0+tshift(:,iatom1,ig)-tshift(:,jatom1,ig)
                      call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     &                            u(2,ig),1)
                      call zone1_number(v1,gbas,ndiv,ir1)
                      jatom0=ip(jatom1,ig)
                      jat0=iat_sort_back(jatom0)
                      j0=nd1*(nat2*(ie1-1)+jat0-1)+1
                      iatom0=ip(iatom1,ig)
                      iat0=iat_sort_back(iatom0)
                      i0=nfun_els_adr(ie,li,isort)+nd*(iat0-1)-nz1+1
                      call dgemm('n','t',nd,nd1,nd1,1.d0,
     &                           g_mm(i0,j0,ir1),nfun1,u(dli1,ig),nd1,
     &                           0.d0,aa,nd)
                      call dgemm('n','n',nd,nd1,nd,1.d0,u(dli,ig),nd,aa,
     &                           nd,1.d0,bb,nd)
                    enddo   !! over ig
                    bb=bb/ngroup
                    do j=1,nd1
                      m1=j-li1
                      lm1=l1*(l1+1)+m1+1
                      km1=indbasa(ii1,ie11,lm1,jsort)
                      do i=1,nd
                        m=i-li
                        lm=l*(l+1)+m+1
                        km=indbasa(ii,ie00,lm,isort)
                        gr(km,km1,ind_ir0)=bb(i,j)
                      enddo
                    enddo
                  enddo  !! over ii
                enddo  !! over ie00
                deallocate(aa,bb)
              enddo   !! over l
            enddo   !! over ii1
          enddo  !! over ie11
        enddo  !! over ind_ir0
        deallocate(g_mm)
      enddo  !! over ind_ls
      deallocate(adr_ls,nd_nls,nm_nls,gx_tau,tmp1,z_ord,trf)
      if(npr/=1) call dgop(gr,nf1*nf2*nd_nrr,'  +',comm_b)
      end
