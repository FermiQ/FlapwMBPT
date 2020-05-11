      subroutine group_heg
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer :: nt,ngen,maxbv,lmm,jmm,ig
	integer, allocatable :: figen(:)
	real*8, allocatable :: fac(:),su1(:,:),gen(:,:),sgen(:,:),ih(:),
     &	                     ihh(:),sig(:)
	complex*16, allocatable :: su1_j(:,:),genj(:,:),ihj(:),ihhj(:),
     &	                       sigj(:)
      inv=0
	allocate(fac(0:8))
      call facc(fac,8)
      nt=index(symgen,' ')-1
      maxbv=1
      maxwig=(maxbv+1)*(2*maxbv+1)*(2*maxbv+3)/3
      maxwig_j=4*(maxbv+1)*(2*(maxbv+1)**2+1)/3
      allocate(u(maxwig,48))
      allocate(uj(maxwig_j,48))
      allocate(ifia(48))
      allocate(shift(3,48))
      lmm=maxbv
	jmm=2*maxbv+1
      allocate(genj(maxwig_j,48))
      allocate(gen(maxwig,48))
      allocate(sgen(3,144))
      allocate(figen(48))
      allocate(ih(maxwig))
      allocate(ihj(maxwig_j))
      allocate(su1(-lmm:lmm,-lmm:lmm))
      allocate(su1_j(-jmm:jmm,-jmm:jmm))
      call decode_generators(symgen,nt,gen,genj,sgen,figen,lmm,jmm,ngen,
     &                       ih,ihj,su1,su1_j,fac)
      allocate(ihh(maxwig))
      allocate(sig(maxwig))
      allocate(ihhj(maxwig_j))
      allocate(sigj(maxwig_j))
      call sgroup(gen,genj,sgen,ngen,u,uj,shift,ngroup,48,lmm,jmm,ih,
     &	          ihj,ihh,ihhj,sig,sigj,su1,su1_j,figen)
      deallocate(genj,gen,sgen,figen,ih,ihj,su1,su1_j,ihh,sig,ihhj,sigj)
      symmorph=.true.
c ----------------------------------------------------------------------
      inv=1
      invers=.true.
      if(maswrk) then
       write(iun,*)'Number of elements in symmetry group is ',ngroup
      endif
	deallocate(fac)
c -------- Cartesian rotation matrices for optimization --------------
      do ig=1,ngroup
        u_opt(1,1,ig)=u(10,ig)
        u_opt(2,1,ig)=u(8,ig)
        u_opt(3,1,ig)=u(9,ig)
        u_opt(1,2,ig)=u(4,ig)
        u_opt(2,2,ig)=u(2,ig)
        u_opt(3,2,ig)=u(3,ig)
        u_opt(1,3,ig)=u(7,ig)
        u_opt(2,3,ig)=u(5,ig)
        u_opt(3,3,ig)=u(6,ig)
      enddo	
      call timel('**** GROUP_HEG finished ************')
      end
