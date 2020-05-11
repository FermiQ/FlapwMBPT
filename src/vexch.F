      subroutine vexch(key,v0,bv0)
	  use atom_mod
      use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
	  integer, intent(in) :: key
	  real*8, intent(out) :: v0(maxmt),bv0(maxmtb)
      integer :: isort,ndimv,ntheta,nphi,nangl,ishang,iang,itheta,iphi,
     &           irad,ispin,l,m,lm,isym,mt,i
	  real*8 :: dphi,zcoord,xcoord,ycoord,stheta,wgtl,phi,excup,excdn,
     &	        vxcup,vxcdn,wv,sm,rotot(2),vxcpot,bpot,bx,by,bz,dum,
     &          vx(2),totm(3),gro(3,3),gmod(3,3),romod(3),rolap(2),
     &          dqdall,sinphi,cosphi
      real*8, allocatable :: zl(:),wl(:),vt(:),bt(:),ylm(:),
     &                       drodr(:,:,:),d2rodr2(:,:,:),dspdr(:,:),
     &                       d2spdr2(:,:),d12yl(:,:),gbar(:)
	  allocate(vt(maxmt))
	  vt=0.d0
	  if(magn==2) then
	    allocate(bt(maxmtb))
	    bt=0.d0
	  endif
	  allocate (ylm(limlpb))
c  **** calculate exchange-correlation by direct expansion ****
      do isort=1,nsort
        ndimv=(lmpb(isort)+1)**2
        ntheta=12*lmpb(isort)+3
        allocate(zl(ntheta))
        allocate(wl(ntheta))
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
        if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
          allocate(gbar(0:nrad(isort)))
          gbar=0.d0
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
c  find contribution to XC potential from "r,theta,phi"
            call exch_corr(rotot(2),rotot(1),ilda,nloc,gro,gmod,romod,
     &                     rolap,EXCUP,EXCDN,VXCUP,VXCDN)
            if(key.eq.1) then
              vx(1)=vxcdn
              vx(2)=vxcup
            else
              vx(1)=excdn
              vx(2)=excup
            endif
            if(key_hybrid==1) then
              if(hybrid_type<=2) then
                call exch_corr(rotot(1),rotot(2),6,nloc,gro,gmod,romod,
     &                         rolap,EXCUP,EXCDN,VXCUP,VXCDN)
              else
                call exch_corr(rotot(1),rotot(2),6,3,gro,gmod,romod,
     &                         rolap,EXCUP,EXCDN,VXCUP,VXCDN)
              endif
              if(key.eq.1) then
                vx(1)=vx(1)-alpha_hbr*vxcdn
                vx(2)=vx(2)-alpha_hbr*vxcup
              else
                vx(1)=vx(1)-alpha_hbr*excdn
                vx(2)=vx(2)-alpha_hbr*excup
              endif
            endif
            if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
              dum=sqrt(romod(3)/(rotot(1)+rotot(2)))
              gbar(irad)=gbar(irad)+wgtl*dum
            endif
            if(magn==1) then
              do ispin=1,nspin
                wv=wgtl*vx(ispin)
                do l=0,lmpb(isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    if(sym(lm,isort))then
                      isym=lmsym(lm,isort)
                      mt=indmt(isym,isort,ispin)
                      vt(mt+irad)=vt(mt+irad)+ylm(lm)*wv
                    endif
                  enddo   !!! over m
                enddo   !!! over l
              enddo   !!! over ispin
            else if(magn.eq.2) then
              vxcpot=0.5d0*(vx(1)+vx(2))
              wv=wgtl*vxcpot
              do l=0,lmpb(isort)
                do m=-l,l
                  lm=l*(l+1)+m+1
                  if(sym(lm,isort))then
                    isym=lmsym(lm,isort)
                    mt=indmt(isym,isort,1)
                    vt(mt+irad)=vt(mt+irad)+ylm(lm)*wv
                  endif
                enddo   !!! over m
              enddo   !!! over l
              bpot=0.5d0*(vx(2)-vx(1))
              if(abs(sm).lt.1.d-10) then
                bx=0.d0
                by=0.d0
                bz=0.d0
              else
                bx=bpot*b_ext(1)
                by=bpot*b_ext(2)
                bz=bpot*b_ext(3)
              endif
              do i=1,3
                if(i.eq.1) wv=wgtl*bx
                if(i.eq.2) wv=wgtl*by
                if(i.eq.3) wv=wgtl*bz
                do l=0,lmpb(isort)
                  do m=-l,l
                    lm=l*(l+1)+m+1
                    if(symb(lm,i,isort))then
                      isym=lmsymb(lm,i,isort)
                      mt=indmtb(isym,isort)
                      bt(mt+irad)=bt(mt+irad)+wv*ylm(lm)
                    endif
                  enddo   !!! over m
                enddo   !!! over l
              enddo   !!! over i
            endif  !!! only for mag = 2
          enddo   !!! over irad
        enddo   !!! over iang
        deallocate(zl,wl)
        if(key==1.and.(hybrid_type==2.or.hybrid_type==4)) then
#ifdef MPI
          call DGOP(gbar,nrad(isort)+1,'  +',MPI_COMM_WORLD)
#endif
          do irad=0,nrad(isort)
            gbar(irad)=gbar(irad)*dr(irad,isort)*r(irad,isort)**2
          enddo
          gbar_param=gbar_param+dqdall(h(isort),gbar,nrad(isort))
     &                         *nhsort(isort)/amega
          deallocate(gbar)
        endif
        if(nloc.gt.0) then
          deallocate(drodr,d2rodr2,d12yl)
          if(magn==2) deallocate(dspdr,d2spdr2)
        endif
      enddo   !!! over isort
#ifdef MPI
      call DGOP(vt,maxmt,'  +',MPI_COMM_WORLD)
      if(magn.eq.2) call DGOP(bt,maxmtb,'  +',MPI_COMM_WORLD)
#endif      
      do i=1,maxmt
        v0(i)=v0(i)+vt(i)
      enddo
      if(magn.eq.2) then
        do i=1,maxmtb
          bv0(i)=bv0(i)+bt(i)
        enddo
      endif
      do isort=1,nsort
        do ispin=1,nspin
          mt=indmt(1,isort,ispin)
          v0(mt)=v0(mt+1)
        enddo
        if(magn.eq.2) then
          do i=1,3
            if(symb(1,i,isort)) then
              isym=lmsymb(1,i,isort)
              mt=indmtb(isym,isort)
              bv0(mt)=bv0(mt+1)
            endif
          enddo
        endif
      enddo   !!! over isort
	  deallocate (ylm)
	  deallocate(vt)
	  if(magn==2) deallocate(bt)
      end
