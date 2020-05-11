      subroutine check_band_states_red(zb,ev)
	  use atom_mod
	  use manager_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: zb(nfun_red,nbndf_bnd,nqdiv_c,nspin),
     &                          ev(nbasmpw_red,nbndf_bnd,nqdiv_c,nspin)
      integer :: igb,isort,iatom,lm,ii,l,m,k,k1,ibnd,n,ipp,ih,is1,
     &           limj,ll,ll0,ist,i,j,mj,lmj,mjs,lmm,i11,lget,li,i_len,
     &           ispin,nba,nbnd,ind,k2,ka,ib,j1,j2
      real*8 :: pi4om,pi2,pi2a,gki(3),gi,sc,cs,sn,fr,sq2,fii,ffr,ee,
     &          ffi,drr,di,bdr,bdi
      complex*16 :: cc,ccd,cc1,ccd1,ccc
	  real*8, allocatable :: ylm(:),bes(:),besd(:)
      complex*16, allocatable :: mtl(:,:,:,:),itl(:,:,:,:),itl0(:,:,:,:)
	  allocate(ylm((maxb+1)**2),bes(0:maxb+1),besd(0:maxb+1))
	  allocate(itl(nrel*(maxb+1)**2,nsort,2,nbndf_bnd))
	  allocate(itl0(nrel*(maxb+1)**2,nsort,2,nbasmpw_red))
	  allocate(mtl(nrel*(maxb+1)**2,nsort,2,nbndf_bnd))
      pi4om=4.d0*pi/sqrt(amega)
      pi2=2.d0*pi
      pi2a=pi2/par
      sq2=sqrt(2.d0)
	  i_len=len_trim(allfile)
	  if(maswrk) open(3,file=allfile(1:i_len)//'BND_STATES_MATCH')
      do ispin=1,nspin
        if(maswrk) write(3,*)' Ispin = ',ispin
        do k=1,npnt_c
          if(maswrk) write(3,*)' K = ',k
          ka=k_a_from_c(k)
          nba=nbask_red(ka)
          nbnd=n_low_bnd(ka,ispin)
          itl0=(0.d0,0.d0)
          do igb=1,nba/nrel
	        ind=indgb_red(igb,ka)
            gki=(pnt_c(:,k)+gbs(:,ind))*pi2a
            gi=sqrt(dot_product(gki,gki))
            if(gi.lt.1.d-4) then
              ylm(1)=1.d0/sqrt(4.d0*pi)
              if(maxb.gt.0) call dcopy((maxb+1)**2-1,0.d0,0,ylm(2),1)
            else
              call sphharm(gki(1)/gi,gki(2)/gi,gki(3)/gi,maxb,ylm)
            endif
            do isort=1,nsort
              iatom=iat_1(isort)
              call bessr(gi,smt(isort),bes,besd,maxb)
              sc=par*dot_product(gki,tau(:,iatom))
              cs=cos(sc)
              sn=sin(sc)
              ccc=dcmplx(cs,sn)
              if(irel<=1) then
                lm=0
                do l=0,maxb
                  cc=pi4om*bes(l)*(0.d0,1.d0)**l
                  ccd=pi4om*besd(l)*(0.d0,1.d0)**l
                  do m=-l,l
                    lm=lm+1
                    cc1=cc*ylm(lm)
                    ccd1=ccd*ylm(lm)
                    itl0(lm,isort,1,igb)=cc1*ccc
                    itl0(lm,isort,2,igb)=ccd1*ccc
                  enddo   !!! over m
                enddo   !!! over l
              else if(irel==2) then
                do ipp=1,2     !! over s
                  ih=(ipp-1)*nba/nrel+igb
                  if(ipp.eq.1) is1=-1
                  if(ipp.eq.2) is1=1
                  limj=0
                  do l=0,maxb
                    ll=l*(l+2)
                    ll0=l*(l+1)+1
                    ist=-1
                    if(l.eq.0) ist=1
                    do i=ist,1,2        !! over i
                      if(i.eq.-1) ii=1
                      if(i.eq.1) ii=2
                      j=l+l+i
                      do mj=-j,j,2        !! =2*mj in (III-61-1)
                        limj=limj+1
                        lmj=ll+(mj+3)/2
                        mjs=(mj-is1)/2     !! =mj-s in (III-61-1)
                        if(iabs(mjs)>l) cycle
                        lm=ll0+mjs
                        lmm=ll0-mjs
                        if(mjs.lt.0) then
                          fr=ckg2(ii,ipp,lmj)*ylm(lmm)/sq2
                          fii=ckg2(ii,ipp,lmj)*ylm(lm)/sq2
                        else if(mjs.gt.0) then
                          i11=(-1)**mjs
                          fr=ckg2(ii,ipp,lmj)*ylm(lm)*i11/sq2
                          fii=-ckg2(ii,ipp,lmj)*ylm(lmm)*i11/sq2
                        else if(mjs.eq.0) then
                          fr=ckg2(ii,ipp,lmj)*ylm(lm)
                          fii=0.d0
                        endif
                        ffr=fr*rr(l)-fii*ri(l)
                        ffi=fr*ri(l)+fii*rr(l)
                        if(rel_interst) then
                        else
                          drr=pi4om*bes(l)*ffr
                          di=pi4om*bes(l)*ffi
                          bdr=pi4om*besd(l)*ffr
                          bdi=pi4om*besd(l)*ffi
                        endif
                        itl0(limj,isort,1,ih)=dcmplx(drr*cs-di*sn,
     &                                               drr*sn+di*cs)
                        itl0(limj,isort,2,ih)=dcmplx(bdr*cs-bdi*sn,
     &                                               bdr*sn+bdi*cs)
                      enddo   !!! over mj
                    enddo   !!! over i
                  enddo   !!! over l
                enddo   !!! over ipp
              endif
            enddo   !!! over isort
          enddo  !!! over igb
          n=nrel*(maxb+1)**2*nsort*2
	      call zgemm('n','n',n,nbnd,nba,(1.d0,0.d0),itl0,n,
     &               ev(1,1,k,ispin),nbasmpw_red,(0.d0,0.d0),itl,n)
	      mtl=(0.d0,0.d0)
	      do isort=1,nsort
	        iatom=iat_1(isort)
	        do lm=1,nrel*(lmb_red(isort)+1)**2
	          if(irel<=1) l=lget(lm)
	          if(irel==2) call getlimj(lm,l,i,mj,li,0)
	          ii=n_e_apw(l,isort)  !! augmentation energy
              j1=indbasa_red(1,ii,lm,isort)
              j2=indbasa_red(2,ii,lm,isort)
              if(j1==0.or.j2==0) cycle
              k1=io_lem_red(iatom)+j1-1
              k2=io_lem_red(iatom)+j2-1
              do ibnd=1,nbnd
                mtl(lm,isort,1,ibnd)=zb(k1,ibnd,k,ispin)
                mtl(lm,isort,2,ibnd)=zb(k2,ibnd,k,ispin)
              enddo
	        enddo
          enddo
          if(maswrk) then
	        do ibnd=1,nbnd
              ib=ind_bands_bnd(ibnd,ka,ispin)
              ee=e_bnd(ib,ka,ispin)-chem_pot
              write(3,*)' BND STATE = ',ibnd,' E - mu = ',
     &                                  ee*13.6058d0,' eV'
	          do isort=1,nsort
                write(3,*)' Isort = ',isort
	            do lm=1,nrel*(maxb+1)**2
                  write(3,'(i5,4(e12.5,1x))')lm,mtl(lm,isort,1,ibnd),
     &	                                        itl(lm,isort,1,ibnd)
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      if(maswrk) close(3)
	  deallocate(mtl,itl,itl0,bes,besd,ylm)
      end
