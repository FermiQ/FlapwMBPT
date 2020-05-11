      subroutine rr_to_gg_fft_red(frr,fgg,npg,npqg,kb,kbq,pt,ptq)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: npg,kb,npqg,kbq
      real*8, intent(in) :: pt(3),ptq(3)
      complex*16, intent(in) :: frr(nr_full_red,nr_full_red)
      complex*16, intent(out) :: fgg(npg,npqg)
	integer :: i,i0,ind
      real*8 :: v(3),pi2,fc
      complex*16 :: cc
      complex*16, allocatable :: tc(:,:),tc1(:,:)
c ----------------------------------------------------------------
      pi2=pi+pi
      allocate(tc(nr_full_red,nr_full_red))
c ----- Multiplication with the phase factor ------------------
      do i=1,nr_full_red
        v(:)=r_pnt_red(1,i)*r0b0_red(:,1)
     &      +r_pnt_red(2,i)*r0b0_red(:,2)
     &      +r_pnt_red(3,i)*r0b0_red(:,3)
        fc=-pi2*dot_product(pt,v)
        cc=exp(dcmplx(0.d0,fc))
        tc(i,:)=frr(i,:)*cc
      enddo
      call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,tc,
     &          -1)
      allocate(tc1(nr_full_red,npg))
      do i=1,npg
        i0=indpw_gw_red(i,kb)
        call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
        tc1(:,i)=tc(ind,:)
      enddo
      tc1=tc1*sqrt(amega)/nr_full_red
      deallocate(tc)
c ----- Multiplication with the phase factor ------------------
      do i=1,nr_full_red
        v(:)=r_pnt_red(1,i)*r0b0_red(:,1)
     &      +r_pnt_red(2,i)*r0b0_red(:,2)
     &      +r_pnt_red(3,i)*r0b0_red(:,3)
        fc=pi2*dot_product(ptq,v)
        cc=exp(dcmplx(0.d0,fc))
        tc1(i,:)=tc1(i,:)*cc
      enddo
      call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),npg,tc1,1)
      do i=1,npqg
        i0=indpw_gw_red(i,kbq)
        call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
        fgg(:,i)=tc1(ind,:)*sqrt(amega)/nr_full_red
      enddo
      deallocate(tc1)
      end