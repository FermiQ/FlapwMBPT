      subroutine fifi_j_prepare(fifi_j,pn,npn)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: npn
      real*8, intent(in) :: pn(3,npn)
      real*8, intent(out) :: fifi_j(maxlfun,maxlfun,0:2*maxb,nspin,
     &                              nspin,nsort,npn)
      integer :: isort,ir,l,mt,mt1,iq,ispin,jspin,lf,lf1,isp,jsp
      real*8 :: pi4,qq,rv,dqdall,q(3),c2,sqpi4,pp,pi2a
      real*8, allocatable :: bess(:,:),bessd(:),work(:),work1(:)
      allocate (bess(0:2*maxb,0:maxnrad),
     &	        bessd(0:2*maxb),work(0:maxnrad),work1(0:maxnrad))
      pi4=4.d0*pi
      pi2a=(pi+pi)/par
      sqpi4=sqrt(pi4)
      c2=clight**2
      do iq=1,npn
        q=pn(:,iq)
        qq=sqrt(dot_product(q,q))*pi2a
        do isort=1,nsort
          do ir=0,nrad(isort)
            rv=r(ir,isort)
            call BESSR(qq,rv,bess(0,ir),bessd,2*lmb(isort))
            rv=rv*rv
            work(ir)=rv*dr(ir,isort)
          enddo
          do l=0,2*lmb(isort)
            do jspin=1,nspin
              jsp=min(jspin,nspin_1)
              do ispin=1,nspin
                isp=min(ispin,nspin_1)
                do lf1=1,lfun(isort)
                  mt1=ind_wf(lf1,isort)
                  do lf=1,lfun(isort)
                    mt=ind_wf(lf,isort)
                    do ir=0,nrad(isort)
                      pp=gfun(mt+ir,isp)*gfun(mt1+ir,jsp)
                      if(irel>=1) pp=pp+gfund(mt+ir,isp)
     &	                             *gfund(mt1+ir,jsp)/c2
                      work1(ir)=work(ir)*bess(l,ir)*pp
                    enddo
                    fifi_j(lf,lf1,l,ispin,jspin,isort,iq)=
     &                     pi4*dqdall(h(isort),work1,nrad(isort))
                  enddo   !! over lf
                enddo   !! over lf1
              enddo   !! over ispin
            enddo    !! over jspin
          enddo  !! over l
        enddo !! over isort
      enddo   !! over iq
      deallocate (bess,bessd,work,work1)
      end
      
      
      
      
      subroutine fifi_j_prep(fifi_j,q,ln,isort,ispin)
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ln,isort,ispin
      real*8, intent(in) :: q(3)
      real*8, intent(out) :: fifi_j(ln,ln,0:2*lmb(isort))
      integer :: ir,l,mt,mt1,lf,lf1
      real*8 :: pi4,qq,rv,dqdall,c2,sqpi4,pp,pi2a
      real*8, allocatable :: bess(:,:),bessd(:),work(:),work1(:)
      allocate (bess(0:2*maxb,0:maxnrad),
     &	      bessd(0:2*maxb),work(0:maxnrad),work1(0:maxnrad))
      pi4=4.d0*pi
      pi2a=(pi+pi)/par
      sqpi4=sqrt(pi4)
      c2=clight**2
      qq=sqrt(dot_product(q,q))*pi2a
      do ir=0,nrad(isort)
        rv=r(ir,isort)
        call BESSR(qq,rv,bess(0,ir),bessd,2*lmb(isort))
        rv=rv*rv
        work(ir)=rv*dr(ir,isort)
      enddo
      do l=0,2*lmb(isort)
        do lf1=1,lfun(isort)
          mt1=ind_wf(lf1,isort)
          do lf=1,lfun(isort)
            mt=ind_wf(lf,isort)
            do ir=0,nrad(isort)
              pp=gfun(mt+ir,ispin)*gfun(mt1+ir,ispin)
              if(irel>=1) pp=pp+gfund(mt+ir,ispin)
     &	                     *gfund(mt1+ir,ispin)/c2
              work1(ir)=work(ir)*bess(l,ir)*pp
            enddo
            fifi_j(lf,lf1,l)=pi4*dqdall(h(isort),work1,nrad(isort))
          enddo   !! over lf
        enddo   !! over lf1
      enddo  !! over l
      deallocate (bess,bessd,work,work1)
      end
