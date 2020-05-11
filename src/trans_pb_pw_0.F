      subroutine trans_pb_pw_0(ka,tran,npb,ind_gg)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ka,npb,ind_gg
      complex*16, intent(out) :: tran(npb)
      integer :: l,iatom,i0,isort,irad,iei,lfi,ind,m,lm,i,k0,
     &           igj,ind_gj,ia,ib,ic,j
      real*8 :: pi2,pi4,pi2a,sqpi4,qg(3),qg2,qg1,fas,rv,dqdall
      complex*16 :: cf,eiph,cf1,cf2
      real*8, allocatable :: ylm(:),work(:),work1(:),work2(:),bess(:,:),
     &                       bessd(:)
      complex*16, allocatable :: sig(:)
      k0=i_kref(ka)
      allocate (ylm(limlpb),sig(0:maxpb),work(0:maxnrad),
     &          work1(0:maxnrad),work2(0:maxnrad),
     &          bess(0:maxpb,0:maxnrad),bessd(0:maxpb))
      pi2=pi+pi
      pi4=pi2+pi2
      pi2a=pi2/par
      sqpi4=sqrt(pi4)
      eiph=(0.d0,1.d0)
      sig(0)=(1.d0,0.d0)
      do l=1,maxpb
        sig(l)=sig(l-1)*eiph
      enddo
      sig=pi4*sig
      qg=pi2a*(pnt(:,ka)+gbs(:,ind_gg))
      qg2=dot_product(qg,qg)
      qg1=sqrt(qg2)
      if(qg1.gt.1.d-6) then
        call sphharm(qg(1)/qg1,qg(2)/qg1,qg(3)/qg1,maxpb,ylm)
      else
        ylm(2:limlpb)=0.d0
        ylm(1)=1.d0/sqpi4
      endif
      do iatom=1,natom
        i0=iopb(iatom)-1
        isort=is(iatom)
        fas=par*dot_product(qg,tau(:,iatom))
        cf=dcmplx(cos(fas),sin(fas))/sqrt(amega)
        do irad=0,nrad(isort)
          rv=r(irad,isort)
          call BESSR(qg1,rv,bess(0,irad),bessd,lmpb(isort))
          rv=rv*rv
          work(irad)=rv*dr(irad,isort)
        enddo
        do l=0,lmpb(isort)
          cf1=cf*sig(l)
          do irad=0,nrad(isort)
            work1(irad)=work(irad)*bess(l,irad)
          enddo
          do iei=1,ntle_pb(l,isort)
            lfi=ind_prod0(iei,l,isort)
            ind=ind_pb(lfi,isort)
            do irad=0,nrad(isort)
              work2(irad)=work1(irad)*fun_pb(ind+irad)
            enddo
            cf2=cf1*dqdall(h(isort),work2,nrad(isort))
            do m=-l,l
              lm=l*(l+1)+m+1
              i=i0+ind_prod(iei,lm,isort)
              tran(i)=cf2*ylm(lm)
            enddo
          enddo
        enddo
      enddo  !! over iatom
      do igj=1,nplwgw(k0)
        i=n_pbmt+igj
        ind_gj=indpw_gw(igj,ka)
        ia=igbs(1,ind_gg)-igbs(1,ind_gj)
        ib=igbs(2,ind_gg)-igbs(2,ind_gj)
        ic=igbs(3,ind_gg)-igbs(3,ind_gj)
        j=indplw(ia,ib,ic)
        if(complex_ro) then
          tran(i)=dcmplx(sovr(j),sovi(j))
        else
          tran(i)=dcmplx(sovr(j),0.d0)
        endif
      enddo
      deallocate(ylm,sig,work,work1,work2,bess,bessd)
      end
