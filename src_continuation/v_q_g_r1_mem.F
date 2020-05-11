      subroutine v_q_g_r1_mem(w_c,khbr)
C Calculates V(q;G;r') where q=1,npnt+par; G=1,nplw_gw; r'=1,nr_full -
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: khbr
      complex*16, intent(out) :: w_c(nplw_gw,ndim3_k(me_k+1),nr_full)
      integer :: k,ind_k,j,j0,i,ind
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),w(:,:)
      allocate(tmp(nr_full,nplw_gw))
      pi2=pi+pi
      w_c=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        allocate(w(nplwgw(k),nplwgw(k)))
        call v_coul_full_ii(nplwgw(k),w,ind_k,khbr)
        tmp=(0.d0,0.d0)
        do j=1,nplwgw(k)
          j0=indpw_gw(j,k)
          call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
          do i=1,nplwgw(k)
            tmp(ind,i)=w(i,j)
          enddo
        enddo
        deallocate(w)
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nplwgw(k),tmp,-1)
c ----- Multiplication with the phase factor ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          fc=-pi2*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
          w_c(1:nplwgw(k),ind_k,i)=tmp(i,1:nplwgw(k))*cc
        enddo
      enddo   !! over ind_k
      deallocate(tmp)
      end
