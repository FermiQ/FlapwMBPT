      subroutine bnd_interp(a,np,kp,n,eb,deb)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: np,n
      real*8, intent(in) :: kp(3,np)
      complex*16, intent(in) :: a(nstar_r,n*nspin)
      real*8, intent(out) :: eb(n,np,nspin),deb(n,np,nspin,3)
      integer :: ispin,k,i,ind,m,j
      real*8 :: pi2,phase
      complex*16 :: cc
      complex*16, allocatable :: sm(:),dsm(:,:)
      pi2=pi+pi
      eb=0.d0
      deb=0.d0
      allocate(sm(2:nstar_r))
      allocate(dsm(3,2:nstar_r))
      do ispin=1,nspin
        do k=1,np
          sm=(0.d0,0.d0)
          dsm=(0.d0,0.d0)
          do m=2,nstar_r
            do i=indstar_r(m-1)+1,indstar_r(m)
              phase=pi2*dot_product(rbs(:,i),kp(:,k))
              cc=dcmplx(cos(phase),sin(phase))
              sm(m)=sm(m)+cc
              cc=cc*(0.d0,1.d0)*par
              dsm(:,m)=dsm(:,m)+cc*rbs(:,i)
            enddo
            sm(m)=sm(m)/(indstar_r(m)-indstar_r(m-1))
            dsm(:,m)=dsm(:,m)/(indstar_r(m)-indstar_r(m-1))
          enddo
          do i=1,n
            ind=(ispin-1)*n+i
            eb(i,k,ispin)=a(1,ind)
            do m=2,nstar_r
              eb(i,k,ispin)=eb(i,k,ispin)+a(m,ind)*sm(m)
              do j=1,3
                deb(i,k,ispin,j)=deb(i,k,ispin,j)+a(m,ind)*dsm(j,m)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(sm,dsm)
      end
