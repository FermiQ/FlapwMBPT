      subroutine g_mm_sort_c(gc,it,ind_tau,ispin,nf1,nf2,nd_nrr,
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
      complex*16, intent(out) :: gc(nf1,nf2,nd_nrr)
      integer :: k0,nbnd,ind,j,ib,ind_k,ig,i_tau,j0,ir,ir0,ieapw,
     &           i0,ind_ir0,i,m,nls,l,npr,ndim0,dli,dli1,
     &           nd1,nd,ind_ls,ls,l1,li1,ne1,ns1,indls1,li,ie,
     &           ie00,ie11,iend,jend,ii,jatom1,jat1,iatom1,
     &           iat1,ie1,ii1,ir1,jatom0,jat0,iatom0,iat0,m1,lm1,
     &           km1,lm,km,iii1,jjj1,iii,jjj,nz1
      real*8 :: pi2,tu,tt(3),phase,v0(3),v(3),v1(3)
      integer, allocatable :: adr_ls(:),nd_nls(:),nm_nls(:,:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp1(:,:),tmp2(:,:),g_mm(:,:,:),
     &                           z_ord(:,:,:),trf(:,:),g_tmp(:,:,:),
     &                           aa(:,:),bb(:,:),uu(:)
      allocate(tmp1(nbndf,nbndf))
      allocate(uu(4*(maxb+1)**2))
      pi2=pi+pi
      i_tau=me_t*ndim3_tau+ind_tau-1
      if(it==2) i_tau=n_tau-i_tau
      nz1=io_lem(iat_1(isort))
c -----------------------------------------------------------------
      allocate(adr_ls(2*lmb(jsort)+1))
      ind=0
      do li=1,2*lmb(jsort)+1
        ind=ind+1
        adr_ls(ind)=li
      enddo
      nls=ind
      npr=nproc_b
      ndim0=nls/npr
      if(ndim0*npr<nls) ndim0=ndim0+1
      allocate(nd_nls(npr))
      allocate(nm_nls(ndim0,npr))
      call size_shift_par_load(nls,npr,nd_nls,nm_nls,ndim0)
      allocate(trf(ndim3_k(me_k+1),nqdiv))
      do ir=1,nqdiv
        tt=rr_3(:,ir)
        do ind_k=1,ndim3_k(me_k+1)
          k0=n3_mpi_k(me_k+1)+ind_k
          phase=pi2*dot_product(pnt(:,k0),tt)
          trf(ind_k,ir)=wgt(k0)*dcmplx(cos(phase),sin(phase))
        enddo
      enddo
c -----------------------------------------------------------------
      gc=(0.d0,0.d0)
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
        li1=adr_ls(ls)
        l1=li1/2
        iii1=2*(li1-l1-l1+2)-5
        jjj1=l1+l1+iii1
        nd1=((li1+1)/2)*2
        ne1=ntle(l1,jsort)
        ieapw=n_e_apw(l1,jsort)
        if(augm(ieapw,l1,jsort)=='LAW'.or.
     &     augm(ieapw,l1,jsort)=='AWL') ne1=ne1+1
        ns1=nd1*ne1*nat2
        indls1=nfun_ls_adr(li1,jsort)
        allocate(g_tmp(nfun1,ns1,ndim3_k(me_k+1)))
        allocate(tmp2(nbndf,ns1))
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
     &               z_ord(nz1,1,ind_k),nfun,tmp2,nbndf,(0.d0,0.d0),
     &               g_tmp(1,1,ind_k),nfun1)
        enddo   !! over ind_k
        deallocate(tmp2)
        allocate(g_mm(nfun1,ns1,nqdiv))
        call zgemm('n','n',nfun1*ns1,nqdiv,ndim3_k(me_k+1),
     &             (1.d0,0.d0),g_tmp,nfun1*ns1,trf,ndim3_k(me_k+1),
     &             (0.d0,0.d0),g_mm,nfun1*ns1)
        deallocate(g_tmp)
        if(nproc_k/=1) call dgop(g_mm,2*nfun1*ns1*nqdiv,'  +',comm_k)
c ------------------------------------------------------------------
        dli1=(8*l1**3+12*l1*l1+4*l1+3)/3+2*l1*l1*(iii1-1)
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
              do li=1,2*lmb(isort)+1
                l=li/2
                iii=2*(li-l-l+2)-5
                jjj=l+l+iii
                nd=((li+1)/2)*2
                allocate(aa(nd,nd1))
                allocate(bb(nd,nd1))
                dli=(8*l**3+12*l*l+4*l+3)/3+2*l*l*(iii-1)
                ie=0
                do ie00=1,ntle(l,isort)
                  iend=1
                  if(augm(ie00,l,isort)=='LAW') iend=2
                  if(augm(ie00,l,isort)=='AWL') iend=2
                  do ii=1,iend
                    ie=ie+1
                    bb=(0.d0,0.d0)
                    do ig=1,ngroup
                      v=v0+tshift(:,iatom1,ig)-tshift(:,jatom1,ig)
                      call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),
     &                            u(2,ig),1)
                      call zone1_number(v1,gbas,ndiv,ir1)
                      jatom0=ip(jatom1,ig)
                      jat0=iat_sort_back(jatom0)
                      j0=nd1*(nhsort(jsort)*(ie1-1)+jat0-1)+1
                      iatom0=ip(iatom1,ig)
                      iat0=iat_sort_back(iatom0)
                      i0=nfun_els_adr(ie,li,isort)+nd*(iat0-1)-nz1+1
                      call zgemm('n','t',nd,nd1,nd1,(1.d0,0.d0),
     &                           g_mm(i0,j0,ir1),nfun1,uj(dli1,ig),nd1,
     &                           (0.d0,0.d0),aa,nd)
                      call zcopy(nd*nd,uj(dli,ig),1,uu,1)
                      uu=conjg(uu)
                      call zgemm('n','n',nd,nd1,nd,(1.d0,0.d0),uu,nd,aa,
     &                           nd,(1.d0,0.d0),bb,nd)
                    enddo   !! over ig
                    bb=bb/ngroup
                    j=0
                    do m1=-jjj1,jjj1,2
                      j=j+1
                      call getlimj(lm1,l1,iii1,m1,li1,1)
                      km1=indbasa(ii1,ie11,lm1,jsort)
                      i=0
                      do m=-jjj,jjj,2
                        i=i+1
                        call getlimj(lm,l,iii,m,li,1)
                        km=indbasa(ii,ie00,lm,isort)
                        gc(km,km1,ind_ir0)=bb(i,j)
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
      deallocate(adr_ls,nd_nls,nm_nls,gx_tau,tmp1,z_ord,trf,uu)
      if(npr/=1) call dgop(gc,2*nf1*nf2*nd_nrr,'  +',comm_b)
      end
