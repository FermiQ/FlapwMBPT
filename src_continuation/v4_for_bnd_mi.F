      subroutine v4_for_bnd_mi(vmi,v2_q,phf)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
	  real*8, intent(out) :: vmi(n_pbmt_red,nr_full_red,nqdiv_c)
	  integer :: ir,iq0,j,i,iq,k11,ka0,ka,np,i0,ind
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:)
c ----------------------------------------------------------------
      vmi=0.d0
      allocate(v2r(nqdiv_c,n_pbmt_red,nr_full_red))
      allocate(tmp(nr_full_red,n_pbmt_red))	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        ka0=k_a_from_c(iq0)
        ka=k_a_from_c(iq)
	    np=nplwgw_red(ka0)
        tmp=(0.d0,0.d0)
        do i=1,np
          i0=indpw_gw_red(i,ka)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          tmp(ind,:)=v2_q(1:n_pbmt_red,n_pbmt_red+i,iq)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &            tmp,-1)
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
        do i=1,nr_full_red
          v2r(k11,:,i)=conjg(phf(i,iq))*tmp(i,:)
        enddo
      enddo   !! over iq
      v2r=v2r/sqrt(amega)
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &          v2r,1)
	  v2r=v2r/nqdiv_c
      do ir=1,nqdiv_c
	    do j=1,nr_full_red
	      do i=1,n_pbmt_red
	        vmi(i,j,ir)=v2r(ir,i,j)
	      enddo
	    enddo
      enddo   !! over ir
      deallocate(v2r,tmp)
      end
