      subroutine bands_gradient(n0,n,k,ispin,zb,ab,grad_bnd)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n0,k,ispin,n
      complex*16, intent(in) :: zb(nfun,n0),ab(nbasmpw,n0)
      complex*16, intent(out) :: grad_bnd(n,n,3)
      integer :: l,k0,nb,i,j,m,lm,iatom,isort,lf,lf1,nf,ind,mt1,mt,l1,
     &           ie1,j1,ie,ir,m1,l1m1,km1,km,jj
      real*8 :: pi2a,v(3),a,b,dqdall,td,tr,tt,cr
      complex*16 :: cc
      real*8, allocatable :: ff(:,:),work(:),uxyz(:,:,:)
      complex*16, allocatable :: s(:,:),s3(:,:),uxyzc(:,:,:),c1(:),
     &                           tmp(:,:,:),s1(:,:)
      pi2a=2.d0*pi/par
      k0=i_kref(k)
      nb=nbask(k0)
c ------- Interstitial part -------------------------------------
      allocate(s(nb,nb))
      call s_overlap(0,s,nb,nb,indgb(1,k))
      allocate(s3(nb,nb))
      allocate(s1(nb,n))
      do j=1,3
        do i=1,nb
          v=pnt(:,k)+gbs(:,indgb(i,k))
          cc=dcmplx(0.d0,pi2a*v(j)/amega)
          s3(:,i)=cc*s(:,i)
        enddo
        call zgemm('n','n',nb,n,nb,(1.d0,0.d0),s3,nb,ab,nbasmpw,
     &             (0.d0,0.d0),s1,nb)
        call zgemm('c','n',n,n,nb,(1.d0,0.d0),ab,nbasmpw,s1,nb,
     &             (0.d0,0.d0),grad_bnd(1,1,j),n)
      enddo
      deallocate(s1,s,s3)
c ------- MT part -------------------------------------
      cc=(0.d0,1.d0)/sqrt(2.d0)
      cr=1.d0/sqrt(2.d0)
      allocate(ff(limlb,6))
      lm=0
      do l=0,maxb
          a=(l+l+1)*(l+l+3)
          b=(l+l-1)*(l+l+1)
        do m=-l,l
          lm=lm+1
          ff(lm,1)=-sqrt(dfloat(l+m+1)*(l+m+2)/a)
          ff(lm,2)=sqrt(dfloat(l-m)*(l-m-1)/b)
          ff(lm,3)=sqrt(dfloat(l-m+1)*(l-m+2)/a)
          ff(lm,4)=-sqrt(dfloat(l+m)*(l+m-1)/b)
          ff(lm,5)=sqrt(dfloat(l-m+1)*(l+m+1)/a)
          ff(lm,6)=sqrt(dfloat(l-m)*(l+m)/b)
        enddo
      enddo
      allocate(tmp(nfun,n,3))
      allocate(work(0:maxnrad))
      do iatom=1,natom
        ind=io_lem(iatom)
        isort=is(iatom)
        nf=lfunm(isort)
        allocate(uxyz(nf,nf,3))
        uxyz=0.d0
        do lf1=1,lfun(isort)
	    mt1=ind_wf(lf1,isort)
          l1=li_val(lf1,isort)-1
          ie1=ie_val(lf1,isort)
          j1=in_val(lf1,isort)
          do lf=1,lfun(isort)
	      mt=ind_wf(lf,isort)
            l=li_val(lf,isort)-1
            ie=ie_val(lf,isort)
            j=in_val(lf,isort)
            do ir=0,nrad(isort)
              work(ir)=gfun(mt+ir,ispin)*gfun1(mt1+ir,ispin)
     &                *dr(ir,isort)*r(ir,isort)**2
            enddo
            td=dqdall(h(isort),work,nrad(isort))
            do ir=0,nrad(isort)
              work(ir)=gfun(mt+ir,ispin)*gfun(mt1+ir,ispin)
     &                *dr(ir,isort)*r(ir,isort)
            enddo
            tr=dqdall(h(isort),work,nrad(isort))
            if(l==l1+1) then
              tt=td-l1*tr
              do m1=-l1,l1
                l1m1=l1*(l1+1)+m1+1
                km1=indbasa(j1,ie1,l1m1,isort)
                do m=-l,l
                  lm=l*(l+1)+m+1
                  km=indbasa(j,ie,lm,isort)
                  if(m==m1+1) uxyz(km,km1,1)=uxyz(km,km1,1)
     &                                      +tt*ff(l1m1,1)
                  if(m==m1-1) uxyz(km,km1,2)=uxyz(km,km1,2)
     &                                      +tt*ff(l1m1,3)
                  if(m==m1) uxyz(km,km1,3)=uxyz(km,km1,3)+tt*ff(l1m1,5)
                enddo
              enddo
            endif
            if(l==l1-1) then
              tt=td+(l1+1)*tr
              do m1=-l1,l1
                l1m1=l1*(l1+1)+m1+1
                km1=indbasa(j1,ie1,l1m1,isort)
                do m=-l,l
                  lm=l*(l+1)+m+1
                  km=indbasa(j,ie,lm,isort)
                  if(m==m1+1) uxyz(km,km1,1)=uxyz(km,km1,1)
     &                                      +tt*ff(l1m1,2)
                  if(m==m1-1) uxyz(km,km1,2)=uxyz(km,km1,2)
     &                                      +tt*ff(l1m1,4)
                  if(m==m1) uxyz(km,km1,3)=uxyz(km,km1,3)+tt*ff(l1m1,6)
                enddo
              enddo
            endif
          enddo   !! over lf
        enddo   !! over lf1
c ----- X+iY; X-iY -------->   X;  Y ------------------------------
        do km1=1,nf
          do km=1,nf
            tt=(uxyz(km,km1,1)+uxyz(km,km1,2))/2
            uxyz(km,km1,2)=(uxyz(km,km1,1)-uxyz(km,km1,2))/2
            uxyz(km,km1,1)=tt
          enddo
        enddo
        allocate(uxyzc(nf,nf,3))
        uxyzc=uxyz
        uxyzc(:,:,2)=uxyzc(:,:,2)/(0.d0,1.d0)
        allocate(c1(nf))
        do i=1,3
          do j=1,nf
            c1=uxyzc(j,:,i)
            do lf=1,lfun(isort)
              l=li_val(lf,isort)-1
              ie=ie_val(lf,isort)
              jj=in_val(lf,isort)
              do m=-l,-1
                lm=l*(l+1)+m+1
                km=indbasa(jj,ie,lm,isort)
                lm=l*(l+1)-m+1
                km1=indbasa(jj,ie,lm,isort)
                uxyzc(j,km,i)=cc*(c1(km)-(-1)**(-m)*c1(km1))
              enddo
              do m=1,l
                lm=l*(l+1)-m+1
                km=indbasa(jj,ie,lm,isort)
                lm=l*(l+1)+m+1
                km1=indbasa(jj,ie,lm,isort)
                uxyzc(j,km1,i)=cr*(c1(km)+(-1)**m*c1(km1))
              enddo
            enddo
          enddo
          do j=1,nf
            c1=uxyzc(:,j,i)
            do lf=1,lfun(isort)
              l=li_val(lf,isort)-1
              ie=ie_val(lf,isort)
              jj=in_val(lf,isort)
              do m=-l,-1
                lm=l*(l+1)+m+1
                km=indbasa(jj,ie,lm,isort)
                lm=l*(l+1)-m+1
                km1=indbasa(jj,ie,lm,isort)
                uxyzc(km,j,i)=cc*(-c1(km)+(-1)**(-m)*c1(km1))
              enddo
              do m=1,l
                lm=l*(l+1)-m+1
                km=indbasa(jj,ie,lm,isort)
                lm=l*(l+1)+m+1
                km1=indbasa(jj,ie,lm,isort)
                uxyzc(km1,j,i)=cr*(c1(km)+(-1)**m*c1(km1))
              enddo
            enddo
          enddo
          call zgemm('n','n',nf,n,nf,(1.d0,0.d0),uxyzc(1,1,i),nf,
     &               zb(ind,1),nfun,(0.d0,0.d0),tmp(ind,1,i),nfun)
        enddo   !! over i
        deallocate(uxyz,uxyzc,c1)
      enddo   !!! over iatom
      deallocate(ff,work)
      do i=1,3
        call zgemm('c','n',n,n,nfun,(1.d0,0.d0),zb,nfun,tmp(1,1,i),nfun,
     &             (1.d0,0.d0),grad_bnd(1,1,i),n)
      enddo
      deallocate(tmp)
      end
