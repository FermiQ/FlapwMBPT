      subroutine w_rr_r_r1(wq,wr)
c	Calculates W(R;r;[r']) where R-all; r=1,nr_full; r'=1,nstar(r1) -
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: wq(nplw_gw,nqdiv)
      real*8, intent(out) :: wr(nqdiv,nr_full)
      integer :: k,i,i0,ind,k0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:),wc(:,:)
      allocate(tmp(nr_full))
      pi2=pi+pi
      allocate(wc(nqdiv,nr_full))
      wc=(0.d0,0.d0)
      do k=1,nqdiv
        k0=i_kref(k)
c ------------------------------------------------------------------
        tmp=(0.d0,0.d0)
        do i=1,nplwgw(k0)
          i0=indpw_gw(i,k)
          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
          tmp(ind)=wq(i,k)
        enddo
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
c ----- Multiplication with the phase factor -----------------------
        call zone1_number(pnt(1,k),rb0,ndiv,ind)
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	  	                       +r_pnt(3,i)*r0b0(:,3)
          fc=pi2*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))
          wc(ind,i)=tmp(i)*cc
        enddo
      enddo   !! over k
      deallocate(tmp)
      call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,wc,1)
      fc=sqrt(amega)*nqdiv
      wr=real(wc)/fc
      deallocate(wc)
      end



      subroutine w_rr_r_r1_x(wq,wr)
c	Calculates W(R;r;[r']) where R-all; r=1,nr_full; r'=1,nstar(r1) -
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: wq(nplw_gw,nqdiv)
      complex*16, intent(out) :: wr(nqdiv,nr_full)
      integer :: k,i,i0,ind,k0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(nr_full))
      pi2=pi+pi
      wr=(0.d0,0.d0)
      do k=1,nqdiv
        k0=i_kref(k)
c ------------------------------------------------------
        tmp=(0.d0,0.d0)
        do i=1,nplwgw(k0) 
          i0=indpw_gw(i,k)
          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
          tmp(ind)=wq(i,k)
        enddo
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,tmp,1)
c ----- Multiplication with the phase factor ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &	  	                       +r_pnt(3,i)*r0b0(:,3)
          fc=pi2*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
          call zone1_number(pnt(1,k),rb0,ndiv,ind)
          wr(ind,i)=tmp(i)*cc
        enddo
      enddo   !! over k
      deallocate(tmp)
      call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,wr,1)
      wr=wr/nqdiv
      end
