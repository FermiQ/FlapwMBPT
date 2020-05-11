      subroutine vxint(key,vintr0,bvintr0)
************************************************************
*  The calculation of the exchange-correlation potential   *
*  in Interstitial Region.                                 *
************************************************************
      use etot_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: key
	  complex*16, intent(out) :: bvintr0(3,nplwro)
	  complex*16, intent(out) :: vintr0(nplwro,nspin)
      logical :: jimag,jreal
	  integer :: ispin,im,ia,i,j,ib
      real*8 :: pi2a,pi4a,const,roup,rodn,excup,excdn,vxcup,vxcdn,ddot,
     &          rh,gro(3,3),rolap(2),amod(3),gmod(3,3)
      complex*16 :: cc,cc1
      real*8, allocatable :: rhotm(:,:),rhoim(:),rhors(:,:),rhois(:,:),
     &                       amag(:,:),gradro(:,:,:),grad2ro(:,:,:,:)
	  complex*16, allocatable :: absmag(:),gbar(:),rho(:,:)
	  allocate(rhotm(mmm,nspmag),rhoim(mmm),rhors(mmm,nspmag),
     &	       rhois(mmm,nspmag),gradro(mmm,3,nspmag),
     &         grad2ro(mmm,3,3,nspmag))
	  if(magn.eq.2) then
        allocate(amag(mmm,3))
	    allocate(absmag(mmm))
	    absmag=(0.d0,0.d0)
      endif
      pi2a=2.d0*pi/par
      pi4a=pi2a**2
      const=dble(nspmag)/2.d0
      allocate(rho(nplwro,nspmag))
      if(nspmag==1) then
        rho(:,1)=rointr(:,1)
      else
        if(irel/=2) then
          rho=rointr
        else
          do i=1,nplwro
            rho(i,1)=dot_product(spintr(:,i),b_ext)
            rho(i,2)=0.5d0*(rointr(i,1)+rho(i,1))
            rho(i,1)=0.5d0*(rointr(i,1)-rho(i,1))
          enddo
        endif
      endif
      jimag=.false.
      if(complex_ro) jimag=.true.
      do ispin=1,nspmag
        call fourier(rhotm(1,ispin),rhoim,indstar,rho(1,ispin),1,igbsf,
     &               gbs,0,0,.true.,jimag)
        rhotm(:,ispin)=rhotm(:,ispin)*const
        if(nloc/=0) then
          do ia=1,3
            call fourier(gradro(1,ia,ispin),rhoim,indstar,rho(1,ispin),
     &                   1,igbsf,gbs,ia,0,.true.,jimag)
            do ib=1,3
              call fourier(grad2ro(1,ib,ia,ispin),rhoim,indstar,
     &                     rho(1,ispin),1,igbsf,gbs,ib,ia,.true.,jimag)
            enddo
          enddo
          call dscal(3*mmm,pi2a*const,gradro(1,1,ispin),1)
          call dscal(9*mmm,pi4a*const,grad2ro(1,1,1,ispin),1)
        endif
      enddo   !!! over ispin
      deallocate(rho)
      if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
        allocate(gbar(mmm))
        gbar=(0.d0,0.d0)
      endif
      do im=1,mmm
        rodn=rhotm(im,1)
        roup=rhotm(im,nspmag)
        if(rodn<0.d0.and.roup>0.d0) rodn=0.99d0*roup
        if(roup<0.d0.and.rodn>0.d0) roup=0.99d0*rodn
        if(nloc.gt.0) then
          rolap(1)=0.d0
          rolap(2)=0.d0
          do i=1,3
            gro(i,1)=gradro(im,i,1)
            gro(i,2)=gradro(im,i,nspmag)
            gro(i,3)=gro(i,1)+gro(i,2)
            rolap(1)=rolap(1)+grad2ro(im,i,i,1)
            rolap(2)=rolap(2)+grad2ro(im,i,i,nspmag)
          enddo
          amod(1)=sqrt(ddot(3,gro(1,1),1,gro(1,1),1))
          amod(2)=sqrt(ddot(3,gro(1,2),1,gro(1,2),1))
          amod(3)=sqrt(ddot(3,gro(1,3),1,gro(1,3),1))
          do i=1,3
            gmod(i,1)=0.d0
            gmod(i,2)=0.d0
            gmod(i,3)=0.d0
            do j=1,3
              gmod(i,1)=gmod(i,1)+gro(j,1)*grad2ro(im,i,j,1)
              gmod(i,2)=gmod(i,2)+gro(j,2)*grad2ro(im,i,j,nspmag)
              gmod(i,3)=gmod(i,3)+gro(j,3)*(grad2ro(im,i,j,1)+
     &                                      grad2ro(im,i,j,nspmag))
            enddo
            if(amod(1).gt.1.d-10)gmod(i,1)=gmod(i,1)/amod(1)
            if(amod(2).gt.1.d-10)gmod(i,2)=gmod(i,2)/amod(2)
            if(amod(3).gt.1.d-10)gmod(i,3)=gmod(i,3)/amod(3)
          enddo
        endif
        call exch_corr(ROUP,RODN,ilda,nloc,gro,gmod,amod,rolap,
     &                 EXCUP,EXCDN,VXCUP,VXCDN)
        if(key.eq.1) then
          rhors(im,1)=vxcdn
          rhors(im,nspmag)=vxcup
        else
          rhors(im,1)=excdn
          rhors(im,nspmag)=excup
        endif
        if(key_hybrid==1) then
          if(hybrid_type<=2) then
            call exch_corr(ROUP,RODN,6,nloc,gro,gmod,amod,rolap,
     &                     EXCUP,EXCDN,VXCUP,VXCDN)
          else
            call exch_corr(ROUP,RODN,6,3,gro,gmod,amod,rolap,
     &                     EXCUP,EXCDN,VXCUP,VXCDN)
          endif
          if(key.eq.1) then
            rhors(im,1)=rhors(im,1)-alpha_hbr*vxcdn
            if(nspmag/=1) rhors(im,nspmag)=rhors(im,nspmag)
     &                                    -alpha_hbr*vxcup
          else
            rhors(im,1)=rhors(im,1)-alpha_hbr*excdn
            if(nspmag/=1) rhors(im,nspmag)=rhors(im,nspmag)
     &                                    -alpha_hbr*excup
          endif
        endif
        rhois(im,1)=0.d0
        if(nspmag/=1) rhois(im,nspmag)=0.d0
        if(magn.eq.2) then
          rh=0.5d0*(rhors(im,1)+rhors(im,2))
          rhors(im,2)=0.5d0*(rhors(im,2)-rhors(im,1))
          rhors(im,1)=rh
        endif
        if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
          rh=amod(3)/(roup+rodn)
          gbar(im)=sqrt(rh)
        endif
      enddo   !!! over im
      if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
        call fft3(mdiv(1),mdiv(2),mdiv(3),1,gbar,-1)
        do i=1,nplwro
          im=igbsf(i)
          cc=gbar(im)/mmm
          if(complex_ro) then
            cc1=dcmplx(sovr(i),sovi(i))
          else
            cc1=dcmplx(sovr(i),0.d0)
          endif
          gbar_param=gbar_param+cc*cc1
        enddo
        deallocate(gbar)
      endif
      if(magn.eq.1) then
        do ispin=1,nspin
          jimag=.false.
          if(complex_ro) jimag=.true.
          call fourier(rhors(1,ispin),rhois(1,ispin),indstar,
     &                 vintr0(1,ispin),0,igbsf,gbs,0,0,.true.,jimag)
        enddo
      else if(magn.eq.2) then
        call fourier(rhors(1,1),rhois(1,1),indstar,vintr0(1,1),0,igbsf,
     &               gbs,0,0,.true.,jimag)
        do im=1,mmm
          do ia=1,3
            amag(im,ia)=rhors(im,2)*b_ext(ia)
          enddo
        enddo
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        do ia=1,3
          call dcopy(mmm,0.d0,0,rhois,1)
          call fourier(amag(1,ia),rhois,indstar,bvintr0(ia,:),0,igbsf,
     &                 gbs,0,0,jreal,jimag)
        enddo
      endif   !!! over magn
	  do ispin=1,nspin
        call symscal(vintr0(1,ispin))
	  enddo
      if(magn.eq.2) then
        jreal=.false.
        jimag=.false.
        if(inv.lt.2) jreal=.true.
        if(inv.ne.1) jimag=.true.
        call symvec(bvintr0,jreal,jimag)
      endif  !!! for magn = 2 only
	  deallocate(rhotm,rhoim,rhors,rhois)
	  if(magn.eq.2) then
        deallocate(amag,absmag)
      endif
      end
