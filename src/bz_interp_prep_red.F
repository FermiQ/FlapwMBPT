      subroutine bz_interp_prep_red
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,i,m,iq,k1,n1
      real*8 :: pi2,phase,c1,c2,xx,xx0
      complex*16 :: cc
      real*8, allocatable :: rhom(:)
      complex*16, allocatable :: a(:,:),sm(:),drk(:,:)
      if(nqdiv_c==1) return
      pi2=pi+pi
      m=n_rr_max
      n1=nqdiv_c-1
c ------------------------------------------------------------------
      allocate(rhom(2:m))
      c1=0.75d0
      c2=0.75d0
      xx0=sqrt(dot_product(rbs(:,2),rbs(:,2)))
      do i=2,m
        xx=sqrt(dot_product(rbs(:,i),rbs(:,i)))
        xx=(xx/xx0)**2
        rhom(i)=(1.d0-c1*xx)**2+c2*xx**3
      enddo
c ------------------------------------------------------------------
      allocate(a(n1,n1))
      a=(0.d0,0.d0)
      allocate(sm(nqdiv_c))
      do i=2,m
        do k=1,nqdiv_c
          phase=pi2*dot_product(pnt_c(:,k),rbs(:,i))
          sm(k)=dcmplx(cos(phase),sin(phase))
        enddo
        do k1=2,nqdiv_c
          cc=conjg(sm(k1)-1.d0)
          do k=2,nqdiv_c
            a(k-1,k1-1)=a(k-1,k1-1)+(sm(k)-1.d0)*cc/rhom(i)
          enddo
        enddo
      enddo
      call invers_h(n1,a,n1)
c ------------------------------------------------------------------
      allocate(drk(2:m,2:nqdiv_c))
      drk=(0.d0,0.d0)
      do i=2,m
        do k=1,nqdiv_c
          phase=pi2*dot_product(pnt_c(:,k),rbs(:,i))
          sm(k)=dcmplx(cos(phase),sin(phase))
        enddo
        do k=2,nqdiv_c
          cc=conjg(sm(k)-1.d0)
          do k1=2,nqdiv_c
            drk(i,k1)=drk(i,k1)+a(k-1,k1-1)*cc/rhom(i)
          enddo
        enddo
      enddo
      deallocate(sm,a,rhom)
c ------------------------------------------------------------------
      allocate(bk_bz_red(2:nqdiv_c))
      bk_bz_red=(0.d0,0.d0)
      do i=2,m
        do k=2,nqdiv_c
          bk_bz_red(k)=bk_bz_red(k)+drk(i,k)
        enddo
      enddo
c ------------------------------------------------------------------
      allocate(ckq_bz(2:nqdiv_c,npnt))
      ckq_bz=(0.d0,0.d0)
      do i=2,m
        do iq=1,npnt
          phase=pi2*dot_product(pnt(:,iq),rbs(:,i))
          cc=dcmplx(cos(phase),sin(phase))
          do k=2,nqdiv_c
            ckq_bz(k,iq)=ckq_bz(k,iq)+drk(i,k)*cc
          enddo
        enddo
      enddo
      deallocate(drk)
      end
