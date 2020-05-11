      subroutine theta_int_xc(key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      logical :: jimag
      integer :: ispin,im,ia,i,j,ib,ii,gbs_number,i1,ind
      real*8 :: pi2a,pi4a,const,roup,rodn,ddot,gro(3,3),rolap(2),
     &          amod(3),gmod(3,3),v(3)
      complex*16 :: cc,ss
      real*8, allocatable :: rhotm(:,:),rhoim(:),rhors(:,:),rhois(:,:),
     &                       gradro(:,:,:),grad2ro(:,:,:,:),thet(:)
      complex*16, allocatable :: fxc(:,:)
      allocate(fxc(mmm,n_ixc_0))
      allocate(thet(n_ixc_0))
      theta_int=(0.d0,0.d0)
      allocate(rhotm(mmm,nspmag),rhoim(mmm),rhors(mmm,nspmag),
     &         rhois(mmm,nspmag),gradro(mmm,3,nspmag),
     &         grad2ro(mmm,3,3,nspmag))
      pi2a=2.d0*pi/par
      pi4a=pi2a**2
      const=dble(nspmag)/2.d0
      jimag=.false.
      if(complex_ro) jimag=.true.
      do ispin=1,nspin
        call fourier(rhotm(1,ispin),rhoim,indstar,rointr(1,ispin),
     &               1,igbsf,gbs,0,0,.true.,jimag)
        if(nloc/=0) then
          do ia=1,3
            call fourier(gradro(1,ia,ispin),rhoim,indstar,
     &                   rointr(1,ispin),1,igbsf,gbs,ia,0,.true.,
     &                   jimag)
            do ib=1,3
              call fourier(grad2ro(1,ib,ia,ispin),rhoim,indstar,
     &                     rointr(1,ispin),1,igbsf,gbs,ib,ia,.true.,
     &                     jimag)
            enddo
          enddo
          call dscal(3*mmm,pi2a,gradro(1,1,ispin),1)
          call dscal(9*mmm,pi4a,grad2ro(1,1,1,ispin),1)
        endif
      enddo   !!! over ispin
      do im=1,mmm
        rodn=rhotm(im,1)*const
        roup=rhotm(im,nspmag)*const
        if(rodn<0.d0.and.roup>0.d0) rodn=0.99d0*roup
        if(roup<0.d0.and.rodn>0.d0) roup=0.99d0*rodn
        if(nloc.gt.0) then
          rolap(1)=0.d0
          rolap(2)=0.d0
          do i=1,3
            gro(i,1)=gradro(im,i,1)
            gro(i,2)=gradro(im,i,nspmag)
            gro(i,3)=(gro(i,1)+gro(i,2))*const
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
            if(amod(3).gt.1.d-10)gmod(i,3)=gmod(i,3)*const/amod(3)
          enddo
        endif
        call theta_loc(ilda,roup,rodn,gro,gmod,amod,rolap,thet,key)
        fxc(im,:)=thet
      enddo   !!! over im
      deallocate(rhotm,rhoim,rhors,rhois,thet)
      do ii=1,n_ixc_0
        call fft3(mdiv(1),mdiv(2),mdiv(3),1,fxc(1,ii),-1)
        do i=1,nplwro
          j=igbsf(i)
          cc=fxc(j,ii)/mmm
          do i1=1,nplwdif_pb
            v=gbs(:,i1)+gbs(:,i)
            ind=gbs_number(v)
            if(ncmpl==1) then
              theta_int(i1,ii)=theta_int(i1,ii)+cc*sovr(ind)
            else
              ss=dcmplx(sovr(ind),sovi(ind))
              theta_int(i1,ii)=theta_int(i1,ii)+cc*ss
            endif
          enddo
        enddo
      enddo
      deallocate(fxc)
      end









      subroutine theta_int_xc0(key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
      use etot_mod
	  use manager_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key
	  integer :: i,im,ispin,j,ii
      real*8 :: const,roup,rodn,gro(3,3),gmod(3,3),romod(3),
     &          rolap(2)
	  real*8, allocatable :: thet(:)
      complex*16, allocatable :: wrk(:,:),rhors(:,:)
	  allocate(rhors(mmm,n_ixc_0))
	  allocate(thet(n_ixc_0))
      theta_int=(0.d0,0.d0)
      const=dble(nspin)/2.d0
      allocate(wrk(mmm,nspin))
      do ispin=1,nspin
        wrk=(0.d0,0.d0)
        do i=1,nplwro
          j=igbsf(i)
          wrk(j,ispin)=rointr(i,ispin)
        enddo
        call fft3(mdiv(1),mdiv(2),mdiv(3),1,wrk(1,ispin),1)
      enddo   !!! over ispin
      gro=0.d0
      gmod=0.d0
      romod=0.d0
      rolap=0.d0
      do im=1,mmm
        rodn=wrk(im,1)*const
        roup=wrk(im,nspin)*const
        call theta_loc(ilda,roup,rodn,gro,gmod,romod,rolap,thet,key)
        rhors(im,:)=thet
      enddo   !!! over im
      do ii=1,n_ixc_0
        call fft3(mdiv(1),mdiv(2),mdiv(3),1,rhors(1,ii),-1)
        do i=1,nplwro
          j=igbsf(i)
          theta_int(i,ii)=rhors(j,ii)/mmm
        enddo
      enddo
	  deallocate(rhors,thet,wrk)
      end
