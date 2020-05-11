      subroutine p_ii_q_mem_1(p_c,pgr)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(in) :: p_c(nqdiv,nr_full)
      complex*16, intent(out) :: pgr(nplw_gw,nqdiv)
      integer :: k,k0,i,i0,ind,k1
      real*8 :: v(3),phase,pi2
      complex*16 :: ex
      complex*16, allocatable :: ff(:)
      pi2=pi+pi
      allocate(ff(nr_full))
      do k=1,nqdiv
        call zone1_number(pnt(1,k),rb0,ndiv,k1)
        k0=i_kref(k)
c ----- Phase factors ------------------
        do i=1,nr_full
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          phase=pi2*dot_product(pnt(:,k),v)
          ex=dcmplx(cos(phase),-sin(phase))
          ff(i)=ex*p_c(k1,i)
        enddo
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,ff,-1)
        do i=1,nplwgw(k0)
          i0=indpw_gw(i,k)
          call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
          pgr(i,k)=ff(ind)
        enddo
      enddo   !! over k
      deallocate(ff)
      pgr=pgr/nr_full*sqrt(amega)
      end
