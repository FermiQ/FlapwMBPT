      subroutine eiglin(ispin,nbas,a,b,ndim,yr,yi,ind_y,ar0,ai0,br0,bi0,
     &	                ar1,ai1,br1,bi1)
	  use atom_mod
	  use manager_mod
	  use solid_mod
      implicit none
	  integer, intent(in) :: ispin
      integer :: ind_y(2*limlb,natom),nbas,ndim,ibas0,iatom,
     &	         isort,ibas,nlct,ilc,jbas,jlc,j,i
      real*8 :: yr(max_bs,nbas),yi(max_bs,nbas),ar0(nbas,nlcloc),
     &          ai0(nbas,nlcloc),br0(nbas,nlcloc),bi0(nbas,nlcloc),
     &          ar1(nbas,nbas),ai1(nbas,nbas),br1(nbas,nbas),
     &          bi1(nbas,nbas)
      complex*16 a(ndim,ndim),b(ndim,ndim),aa
	  complex*16, allocatable :: ztmp(:,:)
c     &&& Adding PW-PW and PW-LO NMT contributions &&&&&&&&&&&&&&&&&
      call eiglin_1(ispin,yr,yi,ind_y,nbas,ar0,ai0,ar1,ai1)
c     &&& LO-LO NMT contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      ibas0=0
      do iatom=1,natom
        isort=is(iatom)
        ibas=ibas0
	    nlct=nlctot(isort)
        do ilc=1,nlct
          ibas=ibas+1
          jbas=ibas0
          do jlc=1,ilc
            jbas=jbas+1
            if(irel.le.1) aa=dcmplx(pv(ilc,jlc,iatom,ispin),0.d0)
            if(irel.eq.2) aa=dcmplx(pv(ilc,jlc,iatom,ispin)
     &                             ,pvj(ilc,jlc,iatom))
            a(ibas,jbas)=a(ibas,jbas)+aa
          enddo   !!! over jlc
        enddo   !!! over ilc
        ibas0=ibas0+nlct
      enddo   !!! over iatom
c     &&& PW-LO full contribution &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      allocate(ztmp(nbas,max(nbas,nlcloc)))
      do j=1,nlcloc
	    do i=1,nbas
	      ztmp(i,j)=dcmplx(ar0(i,j),ai0(i,j))
	    enddo
	  enddo
      call zlacpy(' ',nbas,nlcloc,ztmp,nbas,a(nlcloc+1,1),ndim)
      do j=1,nlcloc
	    do i=1,nbas
	      ztmp(i,j)=dcmplx(br0(i,j),bi0(i,j))
	    enddo
	  enddo
      call zlacpy(' ',nbas,nlcloc,ztmp,nbas,b(nlcloc+1,1),ndim)
c     &&& PW-PW full contribution &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      do j=1,nbas
	    do i=1,nbas
	      ztmp(i,j)=dcmplx(ar1(i,j),ai1(i,j))
	    enddo
	  enddo
      call zlacpy(' ',nbas,nbas,ztmp,nbas,a(nlcloc+1,nlcloc+1),ndim)
      do j=1,nbas
	    do i=1,nbas
	      ztmp(i,j)=dcmplx(br1(i,j),bi1(i,j))
	    enddo
	  enddo
      call zlacpy(' ',nbas,nbas,ztmp,nbas,b(nlcloc+1,nlcloc+1),ndim)
	  deallocate(ztmp)
      end
      
      subroutine eiglin_1(ispin,yr,yi,ind_y,nbas,ar0,ai0,ar1,ai1)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
	  integer, intent(in) :: ispin
      integer :: ind_y(2*limlb,natom),iatom,isort,n,nt1,ind,ndimb,lf,
     &	         jbas,nt,nbas,i,ie,lm
      real*8 :: yr(max_bs,nbas),yi(max_bs,nbas),ar0(nbas,nlcloc),
     &          ai0(nbas,nlcloc),ar1(nbas,nbas),ai1(nbas,nbas)
	  real*8, allocatable :: zrl(:,:),ry(:,:),iy(:,:)
c     &&& Adding PW-PW NMT contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        n=lfunm(isort)-nlctot(isort)
        nt=nlctot(isort)
        nt1=nt+1
        allocate(ry(n,nbas))
        allocate(iy(n,nbas))
        do i=1,n
          lm=lm_isz(nt+i,isort)
          lf=lf_isz(nt+i,isort)
          ie=in_val(lf,isort)
          if(ie==1) ind=ind_y(lm,iatom)
          if(ie==2) ind=ind_y(ndimb+lm,iatom)
          ry(i,:)=yr(ind,:)
          iy(i,:)=yi(ind,:)
        enddo
        allocate(zrl(nbas,n))
c     **************   The REAL  part ******************************
        call dgemm('t','n',nbas,n,n,1.d0,ry,n,pv(nt1,nt1,iatom,ispin),
     &             maxel,0.d0,zrl,nbas)
	    if(irel.eq.2) call dgemm('t','n',nbas,n,n,1.d0,iy,n,
     &                           pvj(nt1,nt1,iatom),maxel,1.d0,zrl,nbas)
        call dgemm('n','n',nbas,nbas,n,1.d0,zrl,nbas,ry,n,1.d0,ar1,nbas)
        call dgemm('n','n',nbas,nbas,n,1.d0,zrl,nbas,iy,n,1.d0,ai1,nbas)
c     ********   The IMAGINE  part *******************************
        call dgemm('t','n',nbas,n,n,-1.d0,iy,n,
     &             pv(nt1,nt1,iatom,ispin),maxel,0.d0,zrl,nbas)
        if(irel.eq.2) call dgemm('t','n',nbas,n,n,1.d0,ry,n,
     &                           pvj(nt1,nt1,iatom),maxel,1.d0,zrl,nbas)
        call dgemm('n','n',nbas,nbas,n,1.d0,zrl,nbas,ry,n,1.d0,ai1,nbas)
        call dgemm('n','n',nbas,nbas,n,-1.d0,zrl,nbas,iy,n,1.d0,ar1,
     &             nbas)
        deallocate(ry,iy,zrl)
      enddo   !!! over iatom
c     &&& Adding PW-LO NMT contributions &&&&&&&&&&&&&&&&&&&&&&&&&&&
      jbas=1
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        n=lfunm(isort)-nlctot(isort)
        nt=nlctot(isort)
        nt1=nt+1
        allocate(ry(n,nbas))
        allocate(iy(n,nbas))
        do i=1,n
          lm=lm_isz(nt+i,isort)
          lf=lf_isz(nt+i,isort)
          ie=in_val(lf,isort)
          if(ie==1) ind=ind_y(lm,iatom)
          if(ie==2) ind=ind_y(ndimb+lm,iatom)
          ry(i,:)=yr(ind,:)
          iy(i,:)=yi(ind,:)
        enddo
        if(nt.ne.0) then
          nt1=nt+1
          call dgemm('t','n',nbas,nt,n,1.d0,ry,n,pv(nt1,1,iatom,ispin),
     &               maxel,1.d0,ar0(1,jbas),nbas)
          call dgemm('t','n',nbas,nt,n,-1.d0,iy,n,pv(nt1,1,iatom,ispin),
     &               maxel,1.d0,ai0(1,jbas),nbas)
	      if(irel.eq.2) then
            call dgemm('t','n',nbas,nt,n,1.d0,iy,n,pvj(nt1,1,iatom),
     &                 maxel,1.d0,ar0(1,jbas),nbas)
            call dgemm('t','n',nbas,nt,n,1.d0,ry,n,pvj(nt1,1,iatom),
     &                 maxel,1.d0,ai0(1,jbas),nbas)
	      endif
          jbas=jbas+nt
        endif
        deallocate(ry,iy)
      enddo   !!! over iatom
      end
