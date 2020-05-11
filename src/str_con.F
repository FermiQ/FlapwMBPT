      subroutine str_con(sc,pnt0,khbr,n1,n2,cd1,cd2)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use solid_mod
      use units_mod
      implicit none
      character*1, intent(in) :: cd1,cd2
      integer, intent(in) :: khbr,n1,n2
      real*8, intent(in) :: pnt0(3)
      complex*16, intent(out) :: sc(n1,n2)
      integer :: nlmi1,nlmi2,lmaxst,nlmst,id,jb,lmaxb,ilb0,ib,ip0,ipa,m,
     &           lmaxa,ila0,lget,i,j,max1,max2
      real*8 :: sig,scalp,dlatm
      complex*16 eiphi
      complex*16, allocatable :: dl(:),b1(:,:),b2(:,:)
      if(cd1=='b') then
        max1=maxpb
        nlmi1=limlpb
      else
        max1=maxs
        nlmi1=limlps
      endif
      if(cd2=='b') then
        max2=maxpb
        nlmi2=limlpb
      else
        max2=maxs
        nlmi2=limlps
      endif
      lmaxst=max1+max2
      nlmst=(lmaxst+1)**2
      allocate(dl(nlmst))
      allocate(b1(nlmi1,nlmi2))
      allocate(b2(nlmi1,nlmi2))
c ------ start loop over connecting vectors -----
      do id=1,ndtau
c ---  Get H(l") for pair id with dtau=r(ib)-r(jb) -----------------
        call lat_sums(pnt0,id,lmaxst,dl,khbr)
c        if(khbr==1) call lat_sum_b(pnt0,dtau(1,id),lmaxst,dl,key,keym)
c ------ Put H(l") into SS(l,l') -----------------------------------
        call str_con_matrix(b1,nlmi1,nlmi2,dl,khbr)
c ----------- SS for pair with -dtau -------------------------------
        do jb=1,nlmi2
          ilb0=lget(jb)
          do ib=1,nlmi1
            ila0=lget(ib)
            sig=1-2*mod(ila0+ilb0,2)
            b2(ib,jb)=sig*dconjg(b1(ib,jb))
          enddo
        enddo
c --- sort into arrays b,bb. jb=center, ib=where expanded ----
        do jb=1,natom
          if(cd2=='b') then
            lmaxb=lmpb(is(jb))
            ilb0=iob2(jb)-1
          else
            lmaxb=lmbs(is(jb))
            ilb0=iobs(jb)-1
          endif
          do ib=1,natom
            ip0=ipdtau(ib,jb) !! if>0 means ib <- jb, otherwise ib -> jb
            ipa=iabs(ip0)
            if(ipa.eq.id) then
              sig=1.d0
              if(ip0.lt.0) sig=-1.d0
              scalp=0.d0
              do m=1,3
                dlatm=tau(m,ib)-tau(m,jb)-sig*dtau(m,ipa) ! lat.vec.ne.0
                scalp=scalp+2.d0*pi*dlatm*pnt0(m)
              enddo
              eiphi=dcmplx(cos(scalp),sin(scalp))  !! phase factor
              if(cd1=='b') then
                lmaxa=lmpb(is(ib))
                ila0=iob2(ib)-1
              else
                lmaxa=lmbs(is(ib))
                ila0=iobs(ib)-1
              endif
c -------------- just multiply with possible phase factor ----------
              if(ip0.gt.0) then
                do j=1,(lmaxb+1)**2
                  do i=1,(lmaxa+1)**2
                    sc(i+ila0,j+ilb0)=eiphi*b1(i,j)
                  enddo
                enddo
              else
                do j=1,(lmaxb+1)**2
                  do i=1,(lmaxa+1)**2
                    sc(i+ila0,j+ilb0)=eiphi*b2(i,j)
                  enddo
                enddo
              endif
            endif
          enddo   !! over jb
        enddo     !! over ib
      enddo   !! over id
      deallocate(dl,b1,b2)
      end
