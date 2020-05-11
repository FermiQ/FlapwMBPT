      subroutine bnd_interp_prep(nn,ebnd,a)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nn
      real*8, intent(in) :: ebnd(nn,npnt,nspin)
      complex*16, intent(out) :: a(nstar_r,nn*nspin)
      integer :: ispin,k,info,i,ind,m,n,ns,nstr
      real*8 :: pi2,phase,c1,c2,xx
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: rhom(:)
      complex*16, allocatable :: smk(:,:),x(:,:),b(:,:),aa(:,:)
      pi2=pi+pi
      nstr=nstar_r
      allocate(smk(2:nstr,npnt))
      smk=(0.d0,0.d0)
      do k=1,npnt
        do m=2,nstr
          do i=indstar_r(m-1)+1,indstar_r(m)
            phase=pi2*dot_product(rbs(:,i),pnt(:,k))
            smk(m,k)=smk(m,k)+dcmplx(cos(phase),sin(phase))
          enddo
          smk(m,k)=smk(m,k)/(indstar_r(m)-indstar_r(m-1))
        enddo
      enddo
c ------------------------------------------------------------------
      allocate(rhom(2:nstr))
      c1=0.75d0
      c2=0.75d0
      do i=2,nstr
        xx=(rcostar_r(i)/rcostar_r(2))**2
        rhom(i)=(1.d0-c1*xx)**2+c2*xx**3
      enddo
c ------------------------------------------------------------------
      do k=1,npnt-1
        smk(:,k)=smk(:,k)-smk(:,npnt)
      enddo
      n=npnt-1
      allocate(x(2:nstr,n))
      do i=2,nstr
        x(i,:)=smk(i,1:n)/rhom(i)
      enddo
      allocate(aa(n,n))
      call zgemm('c','n',n,n,nstr-1,(1.d0,0.d0),x(2,1),nstr-1,
     &           smk(2,1),nstr-1,(0.d0,0.d0),aa,n)
      deallocate(x)
c ------------------------------------------------------------------
      allocate(b(n,nn*nspin))
      do ispin=1,nspin
        do i=1,nn
          ind=(ispin-1)*nn+i
          do k=1,n
            b(k,ind)=ebnd(i,k,ispin)-ebnd(i,npnt,ispin)
          enddo
        enddo
      enddo
      allocate(ipiv(n))
      call zgesv(n,nn*nspin,aa,n,ipiv,b,n,info)
      deallocate(ipiv,aa)
      ns=nstr-1
      call zgemm('n','n',ns,nn*nspin,n,(1.d0,0.d0),smk(2,1),ns,b,n,
     &           (0.d0,0.d0),a(2,1),nstr)
      deallocate(b)
      do i=2,nstr
        a(i,:)=conjg(a(i,:))/rhom(i)
      enddo
      deallocate(rhom)
      do ispin=1,nspin
        do i=1,nn
          ind=(ispin-1)*nn+i
          a(1,ind)=ebnd(i,npnt,ispin)
          do m=2,nstr
            a(1,ind)=a(1,ind)-a(m,ind)*smk(m,npnt)
          enddo
        enddo
      enddo
      deallocate(smk)
      end
