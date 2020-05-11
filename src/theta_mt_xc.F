      subroutine theta_mt_xc(key)
c     key = 0 - result is provided in bosonic 4-vector
c     key = 1 - the result is provided in spin-spin 2*2 form
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer, intent(in) :: key
      integer :: isort,ndimv,ntheta,nphi,nangl,ishang,iang,itheta,iphi,
     &           irad,i,l,m,lm,isym,mt,nn
      real*8 :: dphi,zcoord,xcoord,ycoord,stheta,wgtl,phi,wv,sm,
     &	        rotot(2),totm(3),gro(3,3),gmod(3,3),romod(3),
     &          rolap(2),dum,sinphi,cosphi
      real*8, allocatable :: zl(:),wl(:),ylm(:),thet(:),drodr(:,:,:),
     &                       d2rodr2(:,:,:),dspdr(:,:),d2spdr2(:,:),
     &                       d12yl(:,:)
      nn=16
      allocate(zl(4*maxpb+nn),wl(4*maxpb+nn))
      theta_mt=0.d0
      allocate (ylm(limlpb))
      allocate(thet(n_ixc_0))
c  **** calculate exchange-correlation by direct expansion ****
      do isort=1,nsort
        ndimv=(lmpb(isort)+1)**2
        ntheta=4*lmpb(isort)+nn
        nphi=ntheta
        nangl=nphi*ntheta
        ishang=nangl/nproc
        if(ishang*nproc.lt.nangl) ishang=ishang+1
        call mklegw(ntheta,zl,wl)
        dphi=2.d0*pi/dble(nphi)
        if(nloc.gt.0) then
          allocate(drodr(0:nrad(isort),nsym(isort),nspin))
          allocate(d2rodr2(0:nrad(isort),nsym(isort),nspin))
          if(magn==2) then
            allocate(dspdr(0:nrad(isort),nsymb(isort)))
            allocate(d2spdr2(0:nrad(isort),nsymb(isort)))
          endif
          if(magn==1) then
            call dro1sub(drodr,d2rodr2,dum,dum,isort)
          else if(magn==2) then
            call dro1sub(drodr,d2rodr2,dspdr,d2spdr2,isort)
          endif
          allocate(d12yl(limlpb,5))
        endif
        do iang=1+ishang*me,min(nangl,ishang*(me+1))
          itheta=iang/nphi
          if(mod(iang,nphi).ne.0) itheta=itheta+1
          iphi=iang-nphi*(itheta-1)
          zcoord=zl(itheta)
          stheta=sqrt(1.d0-zcoord*zcoord)
          wgtl=dphi*wl(itheta)
          phi=dble(iphi-1)*dphi
          sinphi=sin(phi)
          cosphi=cos(phi)
          xcoord=cosphi*stheta
          ycoord=sinphi*stheta
          call sphharm(xcoord,ycoord,zcoord,lmpb(isort),ylm)
          if(nloc.gt.0) call d12sylm(stheta,zcoord,sinphi,cosphi,
     &                               d12yl,lmpb(isort),ylm)
          do irad=1,nrad(isort)
c  Set up total density at the point "r,theta,phi"
            if(nloc==0) then
              call inp_lda(irad,ylm,rotot,isort,totm,sm)
            else
              if(magn==2) then
                call inp_gga(stheta,zcoord,sinphi,cosphi,irad,ylm,d12yl,
     &                       rotot,romod,gro,gmod,rolap,drodr,d2rodr2,
     &                       dspdr,d2spdr2,isort,totm,sm)
              else
                call inp_gga(stheta,zcoord,sinphi,cosphi,irad,ylm,d12yl,
     &                       rotot,romod,gro,gmod,rolap,drodr,d2rodr2,
     &                       dum,dum,isort,totm,sm)
              endif
            endif
c  find contribution to exchange-correlation potential from "r,theta,phi"
            call theta_loc(ilda,rotot(2),rotot(1),gro,gmod,romod,rolap,
     &                     thet,key)
            do i=1,n_ixc_0
              wv=wgtl*thet(i)
              do l=0,lmpb(isort)
                do m=-l,l
                  lm=l*(l+1)+m+1
                  if(sym(lm,isort))then
                    isym=lmsym(lm,isort)
                    mt=indmt(isym,isort,1)
                    theta_mt(mt+irad,i)=theta_mt(mt+irad,i)+ylm(lm)*wv
                  endif
                enddo   !!! over m
              enddo   !!! over l
            enddo   !!! over ispin
          enddo   !!! over irad
        enddo   !!! over iang
        if(nloc.gt.0) then
          deallocate(drodr,d2rodr2,d12yl)
          if(magn==2) deallocate(dspdr,d2spdr2)
        endif
      enddo   !!! over isort
#ifdef MPI
      call DGOP(theta_mt,n_ixc_0*maxmt/nspin,'  +',MPI_COMM_WORLD)
#endif      
      do isort=1,nsort
        do i=1,n_ixc_0
          mt=indmt(1,isort,1)
          theta_mt(mt,i)=theta_mt(mt+1,i)
        enddo
      enddo   !!! over isort
      deallocate (ylm,zl,wl,thet)
      end
