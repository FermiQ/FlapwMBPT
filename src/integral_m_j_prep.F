      subroutine integral_m_j_prep(n,lma,gg,rad,iik)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,lma
      real*8, intent(in) :: gg,rad(0:n)
      real*8, intent(out) :: iik(n,0:3,0:lma)
      integer :: i,p,n_gauss,nt,ii,k,l
      real*8 :: pi2,ha,b,a,tg,wg,r1,bb
      real*8, allocatable :: t_gauss(:),w_gauss(:),bess(:),bessd(:),
     &                       cn(:,:)
      pi2=pi+pi
	  allocate (bess(0:lma),bessd(0:lma),cn(0:3,0:lma))
c ----------- For Gauss integrations on interval [0:h(isort)] ----	
      n_gauss=24
      allocate(t_gauss(n_gauss))
      allocate(w_gauss(n_gauss))
	  call mklegw(n_gauss,t_gauss,w_gauss)
      iik=0.d0
      do i=1,n
        bb=rad(i)-rad(i-1)
        nt=gg*bb/pi2+1
        ha=1.d0/nt
        cn=0.d0
        if(gg<1.d-12) then
          cn(0,0)=1.d0
          cn(1,0)=0.5d0
          cn(2,0)=1.d0/3.d0
          cn(3,0)=0.25d0
        else
          do ii=1,nt
            b=ii*ha
            a=b-ha
            do k=1,n_gauss
              tg=(a+b)/2+(b-a)/2*t_gauss(k)
              wg=(b-a)/2*w_gauss(k)
              r1=rad(i-1)+bb*tg
              call BESSR(gg,r1,bess,bessd,lma)
              do p=0,3
                cn(p,:)=cn(p,:)+wg*bess*tg**p
              enddo
            enddo
          enddo
        endif
        do l=0,lma
          do p=0,3
            iik(i,p,l)=iik(i,p,l)+bb*cn(p,l)
          enddo
        enddo
      enddo
      deallocate(t_gauss,w_gauss,bess,bessd,cn)
      end
      
      
      
      subroutine integral_mj(n,lma,gg,rad,sm,ik)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,lma
      real*8, intent(in) :: gg,rad(0:n),sm
      real*8, intent(out) :: ik(0:lma)
      integer :: i,n_gauss,nt,ii,k,l
      real*8 :: pi2,ha,b,a,tg,wg,r1,bb,aa,rv,ff
      real*8, allocatable :: t_gauss(:),w_gauss(:),bess(:),bessd(:)
      pi2=pi+pi
	  allocate (bess(0:lma),bessd(0:lma))
c ----------- For Gauss integrations on interval [0:h(isort)] ----	
      n_gauss=24
      allocate(t_gauss(n_gauss))
      allocate(w_gauss(n_gauss))
	  call mklegw(n_gauss,t_gauss,w_gauss)
      ik=0.d0
      do i=1,n
        bb=rad(i)-rad(i-1)
        nt=gg*bb/pi2+1
        bb=bb/2
        aa=(rad(i)+rad(i-1))/2
        ha=2.d0/nt
        do ii=1,nt
          b=-1.d0+ii*ha
          a=b-ha
          do k=1,n_gauss
            tg=(a+b)/2+(b-a)/2*t_gauss(k)
            wg=(b-a)/2*w_gauss(k)
            r1=aa+bb*tg
            call BESSR(gg,r1,bess,bessd,lma)
            rv=1-(r1/sm)**2
            do l=0,lma
              ff=r1**(l+2)*rv**nfgtl
              ik(l)=ik(l)+bb*wg*bess(l)*ff
            enddo
          enddo
        enddo
      enddo
      deallocate(t_gauss,w_gauss,bess,bessd)
      end
