      subroutine c_alpha_mt(nt,x,isort)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nt,isort
      complex*16, intent(out) :: x(nt,nt,3)
      integer :: lf1,mt1,lf,mt,irad,km1,lm1,l1,i1,mj1,li1,km,nf,
     &           lm,l,i,mj,li
      real*8 :: dqdall,spx,spy,spz,spxi,spyi,spzi
      real*8, allocatable :: f1(:,:),f2(:,:),work(:)
      nf=lfun(isort)
      allocate(f1(nf,nf))
      allocate(f2(nf,nf))
      allocate(work(0:nrad(isort)))
      do lf1=1,nf
        mt1=ind_wf(lf1,isort)
        do lf=1,nf
          mt=ind_wf(lf,isort)
          do irad=0,nrad(isort)
            work(irad)=gfun(mt+irad,1)*gfund(mt1+irad,1)*dr(irad,isort)
     &                                *r(irad,isort)**2
          enddo
          f1(lf,lf1)=dqdall(h(isort),work,nrad(isort))
          do irad=0,nrad(isort)
            work(irad)=gfund(mt+irad,1)*gfun(mt1+irad,1)
     &                                 *dr(irad,isort)*r(irad,isort)**2
          enddo
          f2(lf,lf1)=dqdall(h(isort),work,nrad(isort))
        enddo  !!! over lf
      enddo  !!! over lf1
      deallocate(work)
      do km1=1,nt
        lm1=lm_isz(km1,isort)
        lf1=lf_isz(km1,isort)
        call getlimj(lm1,l1,i1,mj1,li1,0)
        do km=1,nt
          lm=lm_isz(km,isort)
          lf=lf_isz(km,isort)
          call getlimj(lm,l,i,mj,li,0)
          call sp_rel(l,i,mj,l1+i1,-i1,mj1,spx,spy,spz,spxi,spyi,spzi)
          x(km,km1,1)=f1(lf,lf1)*dcmplx(spx,spxi)
          x(km,km1,2)=f1(lf,lf1)*dcmplx(spy,spyi)
          x(km,km1,3)=f1(lf,lf1)*dcmplx(spz,spzi)
          call sp_rel(l+i,-i,mj,l1,i1,mj1,spx,spy,spz,spxi,spyi,spzi)
          x(km,km1,1)=x(km,km1,1)-f2(lf,lf1)*dcmplx(spx,spxi)
          x(km,km1,2)=x(km,km1,1)-f2(lf,lf1)*dcmplx(spy,spyi)
          x(km,km1,3)=x(km,km1,1)-f2(lf,lf1)*dcmplx(spz,spzi)
        enddo  !!! over km
      enddo  !!! over km1
      deallocate(f1,f2)
      x=(0.d0,2.d0)*x
      end
