      subroutine theta_q_pw(iq,npw,thet_q)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: npw,iq
      complex*16, intent(out) :: thet_q(npw,npw,n_ixc_0)
      integer :: iatom,isort,j,i,j0,i0,gbs_number,ii,ix,iii,ij,ir,
     &           isym,k,l,m,lm,mt
      real*8 :: v(3),v1(3),pi2a,va,fmt,pi4om,rv,phase,sqpi4,dqdall
      complex*16 :: cc
      integer, allocatable :: atas(:,:),indd(:)
      real*8, allocatable :: ylm(:),bessd(:),bess(:,:),work(:),work1(:)
      complex*16, allocatable :: ff(:)
      pi2a=(pi+pi)/par
      pi4om=4.d0*pi/amega
      sqpi4=sqrt(4.d0*pi)
      ii=0
      allocate(atas(npw,npw))
      allocate(indd(npw*npw))
      do j=1,npw
        j0=indpw_gw(j,iq)
        do i=1,j
          i0=indpw_gw(i,iq)
          v=gbs(:,j0)-gbs(:,i0)
          ij=gbs_number(v)
          iii=0
          do k=1,ii
            if(indd(k)==ij) then
              atas(i,j)=k
              iii=1
            endif
          enddo
          if(iii==0) then
            ii=ii+1
            indd(ii)=ij
            atas(i,j)=ii
          endif
        enddo
      enddo
c ------------------------------------------------------------------
      allocate(ff(ii))
      allocate(ylm(limlpb))
      do ix=1,n_ixc_0
        ff=(0.d0,0.d0)
        do i0=1,ii
          i=indd(i0)   !! in GBS list
          v=gbs(:,i)*pi2a
          va=sqrt(dot_product(v,v))
c --------- MT part ------------------------------------------------
          if(va.gt.1.d-6) then
            call sphharm(v(1)/va,v(2)/va,v(3)/va,maxpb,ylm)
          else
            ylm(2:limlpb)=0.d0
            ylm(1)=1.d0/sqpi4
          endif
          allocate(bess(0:maxpb,0:maxnrad))
          allocate(bessd(0:maxpb))
          allocate(work(0:maxnrad))
          allocate(work1(0:maxnrad))
          do iatom=1,natom
            isort=is(iatom)
            do ir=0,nrad(isort)
              rv=r(ir,isort)
              call BESSR(va,rv,bess(0,ir),bessd,lmpb(isort))
              rv=rv*rv
              work(ir)=rv*dr(ir,isort)
            enddo
            do l=0,lmpb(isort)
              do m=-l,l
                lm=l*(l+1)+m+1
                if(sym(lm,isort))then
                  isym=lmsym(lm,isort)
                  mt=indmt(isym,isort,1)
                  do ir=0,nrad(isort)
                    work1(ir)=work(ir)*bess(l,ir)*theta_mt(mt+ir,ix)
                  enddo
                  fmt=dqdall(h(isort),work1,nrad(isort))
                  phase=dot_product(v,tau(:,iatom))*par
                  cc=(0.d0,1.d0)**l*dcmplx(cos(phase),sin(phase))
     &                             *ylm(lm)
                  ff(i0)=ff(i0)+pi4om*cc*fmt
                endif
              enddo  !! over m
            enddo   !! over l
          enddo   !! over iatom
          deallocate(bess,bessd,work,work1)
c --------- Interstitial -------------------------------------------
          do j=1,nplwro
            v1=gbs(:,i)+gbs(:,j)  !! G'-G+G"
            j0=gbs_number(v1)   !! G_i-G_j in GBS-list
            if(complex_ro) then
              ff(i0)=ff(i0)+theta_int(j,ix)*dcmplx(sovr(j0),sovi(j0))
            else
              ff(i0)=ff(i0)+theta_int(j,ix)*dcmplx(sovr(j0),0.d0)
            endif
          enddo
        enddo   !! over i0
        do j=1,npw
          do i=1,j
            ij=atas(i,j)
            thet_q(i,j,ix)=ff(ij)
          enddo
        enddo
        do j=1,npw-1
          do i=j+1,npw
            thet_q(i,j,ix)=conjg(thet_q(j,i,ix))
          enddo
        enddo
      enddo  !! over ix
      deallocate(atas,ylm,ff,indd)
      end
