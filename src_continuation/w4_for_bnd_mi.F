      subroutine w4_for_bnd_mi(vmi,wmi,w_tau,v2_q,phf)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  real*8, intent(in) :: w_tau(n_pbtot_red,n_pbtot_red,
     &                            ndim3_tau,nqdiv_c)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
	  real*8, intent(out) :: vmi(n_pbmt_red,nr_full_red,nqdiv_c),
     &                           wmi(n_pbmt_red,nr_full_red,nqdiv_c,
     &                           0:n_tau/2)
	  integer :: ind_tau,ir,iq0,j,i,i_tau,iq,k11,ii,ka0,ka,np,i0,ind
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:),tt(:,:)
c ----------------------------------------------------------------
      wmi=0.d0
      vmi=0.d0
      allocate(v2r(nqdiv_c,n_pbmt_red,nr_full_red))
      allocate(tt(n_pbtot_red,n_pbtot_red))
      allocate(tmp(nr_full_red,n_pbmt_red))
c ---- Static part ------------------------------------------------	 
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
c ---- Dynamic part ------------------------------------------------
      do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
c --------------------------------------------------------------- 
        v2r=(0.d0,0.d0)
        do iq=1,nqdiv_c
          call boson_unpack_tau(tt,n_pbtot_red,
     &                          w_tau(1,1,ind_tau,iq),n_pbtot_red,
     &                          n_pbtot_red)
          iq0=i_kref_c(iq)
          ka0=k_a_from_c(iq0)
          ka=k_a_from_c(iq)
	      np=nplwgw_red(ka0)
          tmp=(0.d0,0.d0)
          do i=1,np
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            tmp(ind,:)=tt(1:n_pbmt_red,n_pbmt_red+i)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              tmp,-1)
	      call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
          do i=1,nr_full_red
            v2r(k11,:,i)=conjg(phf(i,iq))*tmp(i,:)
          enddo
        enddo   !! over iq
        v2r=v2r/sqrt(amega)
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &            v2r,1)
	    v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
	      do j=1,nr_full_red
	        do i=1,n_pbmt_red
	          wmi(i,j,ir,i_tau)=v2r(ir,i,j)
	        enddo
	      enddo
        enddo   !! over ir
      enddo  !! over ind_tau
      deallocate(tt,v2r,tmp)
	  if(nproc_t/=1) then
	    ii=n_pbmt_red*nr_full_red*nqdiv_c*(n_tau/2+1)
	    call dgop(wmi,ii,'  +',comm_t)
	  endif
      end
