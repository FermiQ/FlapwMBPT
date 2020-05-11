      subroutine lda_k_point(ind_k,ispin,kp,err_val,err_slo,ebn,nbn,zbn,
     &                       evbn,indg,nbas,key)
c     key = 0 - regular use
c     key = 1 - for band plotting
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,indg(nbndtm),nbas,key
      real*8, intent(in) :: kp(3)
      real*8, intent(inout) :: err_val,err_slo
      integer, intent(out) :: nbn
      real*8, intent(out) :: ebn(nbndf)
      complex*16, intent(out) :: zbn(nfun,nbndf),evbn(nbasmpw,nbndf)
      integer :: nev,neig,ibnd,ind1,iatom,isort,nlct,ndimb,i0,ie,ndim,j,
     &           i,iv,k,ind,lm,ii,n,nn,nt,lf,in
      real*8 :: eps_degen,vabs,v(3),pi2a,acm1,acm2,pi2,phase
      integer, allocatable :: ind_y(:,:)
      real*8, allocatable :: yr(:,:),yi(:,:),ar0(:,:),diag(:),
     &                       ar1(:,:),br0(:,:),br1(:,:),ai0(:,:),
     &                       ai1(:,:),bi0(:,:),bi1(:,:),gb(:,:),ac(:),
     &                       amsm(:,:,:)
      complex*16, allocatable :: a(:,:),b(:,:),yc(:,:),ztmp(:,:),
     &                           z_tmp(:),eikt(:)
      pi2=pi+pi
      k=n3_mpi_k(me_k+1)+ind_k
      eps_degen=1.d-7
      allocate(eikt(natom))
      do iatom=1,natom
        phase=pi2*dot_product(kp,tau(:,iatom))
        eikt(iatom)=dcmplx(cos(phase),sin(phase))
      enddo
      allocate(gb(3,nplwbas))
      do i=1,nbas/nrel
        gb(:,i)=gbs(:,indg(i))
      enddo
      ndim=nbas
      pi2a=(pi+pi)/par
      ndim=nbas+nlcloc
      allocate(diag(ndim))
      allocate(a(ndim,ndim),b(ndim,ndim))
      a=(0.d0,0.d0)
      b=(0.d0,0.d0)
      allocate(ind_y(2*limlb,natom))
      allocate(yr(max_bs,nbas))
      allocate(yi(max_bs,nbas))
      yr=0.d0
      yi=0.d0
      call setk2(kp,gb,yr,yi,nbas,ind_y)
      allocate(ar1(nbas,nbas),br1(nbas,nbas),ai1(nbas,nbas),
     &         bi1(nbas,nbas))
c ------ Interstitial H and O --------------------------------------
      call setwmt(gb,kp,nbas,ar1,ai1,br1,bi1,indg,ispin,1)
      if(k.eq.1) call timel('****** Interstitial (k=1) finished *')
      allocate(ar0(nbas,nlcloc),br0(nbas,nlcloc),ai0(nbas,nlcloc),
     &         bi0(nbas,nlcloc))
c     &&& MT spherical contribution H and O &&&&&&&&&&&&&&&&&&&&&&&&
      call set_h0(b,a,nbas,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
     &            ar1,ai1,br1,bi1,ispin)
      if(k.eq.1) call timel('****** MT-spherical (k=1) finished *')
c  -------- MT non-spherical contribution H ------------------------
      call eiglin(ispin,nbas,a,b,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
     &            ar1,ai1,br1,bi1)
      deallocate(ar0,ar1,ai0,ai1,br0,br1,bi0,bi1)
c ----- Correction for the phase factor e^(ikt) --------------------
      i0=0
      do iatom=1,natom
        isort=is(iatom)
        do i=1,nlctot(isort)
          a(nlcloc+1:ndim,i0+i)=a(nlcloc+1:ndim,i0+i)*eikt(iatom)
          b(nlcloc+1:ndim,i0+i)=b(nlcloc+1:ndim,i0+i)*eikt(iatom)
        enddo
        i0=i0+nlctot(isort)
      enddo
c ------------------------------------------------------------------      
      do j=1,ndim
        do i=1,j-1
          a(i,j)=conjg(a(j,i))
          b(i,j)=conjg(b(j,i))
        enddo
      enddo
      if(k.eq.1) call timel('****** MT-FP (k=1) finished ********')
      allocate(ztmp(ndim,ndim))
c -------- Spectrum of complex matrices ---------------------------
      if(me_t_b==0) call eig_val_gener_zhegvx(ndim,ndim,a,b,nev,diag)
      if(k.eq.1) call timel('****** Diag-n (k=1) finished *******')
      if(nproc_t*nproc_b/=1) then
        call brdcst(msgdbl,a,16*ndim**2,0,comm_t_b)
        call brdcst(msgdbl,diag,8*ndim,0,comm_t_b)
        call brdcst(msgint,nev,4,0,comm_t_b)
      endif
      ztmp=a
      neig=nev
      deallocate(a,b)
c     &&&&&& Getting New N_BND(k) &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(nbndf>=neig) then
        nbn=neig
      else if(nbndf<neig) then
        do ibnd=nbndf,1,-1
          if(diag(ibnd+1)-diag(ibnd).gt.1.d-5)then
            nbn=ibnd
            exit
          endif
        enddo
      endif
      neig=nbn      
      nbn=neig
      ebn(1:nbn)=diag(1:nbn)
      if(nbndf>nbn) ebn(nbn+1:nbndf)=ebn(nbn)+0.2d0
      deallocate(diag)
c ------ Gradients with respect to momentum ------------------------
      if(l_grad_k.and.key==0) then
        call de_dk_num
      endif
c -- Copying first NLCLOC components of EV_LDA into Z for LO - basis
      call zlacpy(' ',nlcloc,neig,ztmp,ndim,zbn,nfun)
c ----- Correction for the phase factor e^(ikt) --------------------
      i0=0
      do iatom=1,natom
        isort=is(iatom)
        do i=1,nlctot(isort)
          zbn(i0+i,1:neig)=zbn(i0+i,1:neig)*eikt(iatom)
        enddo
        i0=i0+nlctot(isort)
      enddo
c  ------- Forming rest components of Z from PW-part of basis ------
      nn=nfun-nlcloc
      allocate(yc(nn,nbas))
      ii=0
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        n=lfunm(isort)-nlctot(isort)
        nt=nlctot(isort)
        do i=1,n
          ii=ii+1
          lm=lm_isz(nt+i,isort)
          lf=lf_isz(nt+i,isort)
          in=in_val(lf,isort)
          if(in==1) ind=ind_y(lm,iatom)
          if(in==2) ind=ind_y(ndimb+lm,iatom)
          yc(ii,:)=dcmplx(yr(ind,:),yi(ind,:))
        enddo
      enddo
      call zgemm('n','n',nn,neig,nbas,(1.d0,0.d0),yc,nn,
     &           ztmp(nlcloc+1,1),ndim,(0.d0,0.d0),zbn(nlcloc+1,1),nfun)
      deallocate(yc)
c - Transformation ZC from FLAPW representation to MB representation -
      allocate(z_tmp(nfun))
      do ibnd=1,neig
        z_tmp=zbn(:,ibnd)
c --------- Local orbitals -------------------------------------
        ind1=1
        do iatom=1,natom
          isort=is(iatom)
          nlct=nlctot(isort)
          call zcopy(nlct,z_tmp(ind1),1,zbn(io_lem(iatom),ibnd),1)
          ind1=ind1+nlct
        enddo
c --------- Augmented orbitals -------------------------------------
        ind1=nlcloc+1
        do iatom=1,natom
          isort=is(iatom)
          ndimb=lfunm(isort)-nlctot(isort)
          nlct=nlctot(isort)
          i0=io_lem(iatom)+nlct
          call zcopy(ndimb,z_tmp(ind1),1,zbn(i0,ibnd),1)
          ind1=ind1+ndimb
        enddo
      enddo
      deallocate(z_tmp)
c     &&&&& Forming EV_LDA &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      call zlacpy(' ',nbas,neig,ztmp(nlcloc+1,1),ndim,evbn,nbasmpw)
c ------------------------------------------------------------------
      allocate(ac(nbn))
      allocate(amsm(2,nbn,nsort))
      ac=0.d0
      amsm=0.d0
      if(irel/=2.or..not.rel_interst) 
     &	call check_band_states(zbn,evbn,kp,gb,nbn,nbas,amsm,ispin)
      call norma_wf(zbn,evbn,gb,nbn,nbas,ac,ispin,indg,kp)
      if(maswrk) then
        if(iprn.ge.1) then
          write(iun,'(a21,i4,3(1x,f8.4))')'Band energies at k = ',k,
     &                                    (kp(i),i=1,3)
          if(iheg) then
            do ie=1,nbas
              iv=indg(ie)
              v=kp+gbs(:,iv)
              vabs=dot_product(v,v)*pi2a**2
              write(iun,*)'e_emp e_found ',ie,vabs,ebn(ie)
            enddo   !!! over ie
          else
            do ie=1,nbn
              acm1=maxval(amsm(1,ie,:))
              acm2=maxval(amsm(2,ie,:))
              write(iun,'(a17,i4,f18.13,3(1x,e12.5))')
     &	        'e d_val d_slp NRM',ie,ebn(ie),acm1,acm2,ac(ie)
            enddo   !!! over ie
          endif
        else if(iprn==0.and.k==1) then
          write(iun,*)'Band energies at k = ',1
          do ie=1,nbn
            acm1=maxval(amsm(1,ie,:))
            acm2=maxval(amsm(2,ie,:))
            write(iun,'(a17,i4,f18.13,3(1x,e12.5))')
     &        'e d_val d_slp NRM',ie,ebn(ie),acm1,acm2,ac(ie)
          enddo   !!! over ie
        endif
      endif
      err_val=max(err_val,maxval(amsm(1,:,:)))
      err_slo=max(err_slo,maxval(amsm(2,:,:)))
      deallocate(ac,amsm) 
      deallocate(ztmp)
      deallocate(gb,yr,yi,ind_y,eikt)
      end
