      subroutine group(lmm,jmm,maxw2)
****************************************************************
*     Program for preparation of transformation matrices       *
*     for spherical harmonics according to the point group.    *
*     To find spherical harmonic which after applying the      *
*     declared group operation will reduce to y(l) multiply    *
*     vector y(m) to matrice u(m,m1).                          *
*      input:                                                  *
*        rbas - translational vectors                          *
*        lmm - lmax i.e 0:s,1:p,2:d,3:f,etc                    *
*        ngroup - number of group elements founded             *
*        u - wigner"s matrices                                 *
*        shift - shifts if group is non-symmorphic one         *
*        inv = 0 no  inversion                                 *
*        inv = 1 yes inversion, ifia=1                         *
*        inv = 2 yes inversion, ifia=-1                        *
*     note that:                                               *
*            u    (l,g) = u(m'm,l,g) = u(iwig,g)               *
*             m'm                                              *
*     where iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m)+l+m'+1      *
****************************************************************
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
	  integer, intent(in) :: lmm,jmm,maxw2
	  integer :: nt,ngen,ir,l,m,m1,iwig,ig
      real*8 :: dshift
	  integer, allocatable :: figen(:)
	  real*8, allocatable :: fac(:),su1(:,:),gen(:,:),sgen(:,:),ih(:),
     &	                     ihh(:),sig(:)
	  complex*16, allocatable :: su1_j(:,:),genj(:,:),ihj(:),ihhj(:),
     &                           sigj(:)
      inv=0
	  allocate(fac(0:8*maxw2))
      call facc(fac,8*maxw2)
      nt=index(symgen,' ')-1
      allocate(genj(maxwig_j,148))
      allocate(gen(maxwig,148))
      allocate(sgen(3,144))
      allocate(figen(148))
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
      call sgroup(gen,genj,sgen,ngen,u,uj,shift,ngroup,148,lmm,jmm,ih,
     &	          ihj,ihh,ihhj,sig,sigj,su1,su1_j,figen)
      deallocate(genj,gen,sgen,figen,ih,ihj,su1,su1_j,ihh,sig,ihhj,sigj)
      symmorph=.true.
      do ir=1,ngroup
        call zone0_cart(shift(1,ir),rbas,shift(1,ir))
        dshift=shift(1,ir)**2+shift(2,ir)**2+shift(3,ir)**2
        if(dshift.gt.1.d-5) symmorph=.false.
        dshift=0.d0
        do l=0,lmm
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if(m1.eq.m) then
                dshift=dshift+abs(u(iwig,ir)-(-1.d0)**l)
              else
                dshift=dshift+abs(u(iwig,ir))
              endif
            enddo
          enddo
        enddo
        if(dshift.lt.1.d-5) then
          inv=1+(1-ifia(ir))/2
          inv_num=ir
        endif
      enddo
c --------------------------------------------------------------------------
      invers=.false.
      if(inv.eq.1) invers=.true.
      if(maswrk) then
       write(iun,*)'Number of elements in symmetry group is ',ngroup
       if(inv.eq.1) write(iun,*)'Inversional operation is presented'
       if(inv.eq.2) write(iun,*)'Inversional operation is not presented'
       if(symmorph) write(iun,*)'The group is symmorphic'
       if(.not.symmorph) write(iun,*)'The group is not symmorphic'
      endif
	  deallocate(fac)
	  real_things=.false.
	  complex_ro=.true.
	  if(invers) then
c	    if(maxval(abs(shift(:,inv_num)))<1.d-6) complex_ro=.false.
	  endif
	  ncmpl=1
	  if(complex_ro) ncmpl=2
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
c ------ Multiplication table -----------------------------------------	
      call group_table
      if(irel==2) then
        allocate(uujj(2,2,2,2,ngroup))
        do ig=1,ngroup
          uujj(1,1,1,1,ig)=conjg(uj(1,ig))*uj(1,ig)
          uujj(2,1,1,1,ig)=conjg(uj(3,ig))*uj(1,ig)
          uujj(1,2,1,1,ig)=conjg(uj(1,ig))*uj(3,ig)
          uujj(2,2,1,1,ig)=conjg(uj(3,ig))*uj(3,ig)
          uujj(1,1,2,1,ig)=conjg(uj(2,ig))*uj(1,ig)
          uujj(2,1,2,1,ig)=conjg(uj(4,ig))*uj(1,ig)
          uujj(1,2,2,1,ig)=conjg(uj(2,ig))*uj(3,ig)
          uujj(2,2,2,1,ig)=conjg(uj(4,ig))*uj(3,ig)
          uujj(1,1,1,2,ig)=conjg(uj(1,ig))*uj(2,ig)
          uujj(2,1,1,2,ig)=conjg(uj(3,ig))*uj(2,ig)
          uujj(1,2,1,2,ig)=conjg(uj(1,ig))*uj(4,ig)
          uujj(2,2,1,2,ig)=conjg(uj(3,ig))*uj(4,ig)
          uujj(1,1,2,2,ig)=conjg(uj(2,ig))*uj(2,ig)
          uujj(2,1,2,2,ig)=conjg(uj(4,ig))*uj(2,ig)
          uujj(1,2,2,2,ig)=conjg(uj(2,ig))*uj(4,ig)
          uujj(2,2,2,2,ig)=conjg(uj(4,ig))*uj(4,ig)
        enddo
      endif
      call timel('**** GROUP finished ****************')					
      end
