      subroutine w4_for_bnd_mm(v4_r,w4_r,w_tau,v2_q,npb)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: npb
	real*8, intent(in) :: w_tau(npb,npb,ndim3_tau,nqdiv_c)
      complex*16, intent(in) :: v2_q(npb,npb,nqdiv_c)
	real*8, intent(out) :: v4_r(n_pbmt_red,n_pbmt_red,nqdiv_c),
     &                       w4_r(n_pbmt_red,n_pbmt_red,nqdiv_c,
     &                            0:n_tau/2)
	integer :: ind_tau,ir,j,i,i_tau,ii,iq,k11
	complex*16, allocatable :: v2r(:,:,:),tmp(:,:)
c ----------------------------------------------------------------
      w4_r=0.d0
      v4_r=0.d0
      allocate(v2r(nqdiv_c,n_pbmt_red,n_pbmt_red))
      allocate(tmp(n_pbmt_red,n_pbmt_red))
c ---- Static part ------------------------------------------------	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
	  call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	  v2r(k11,:,:)=v2_q(1:n_pbmt_red,1:n_pbmt_red,iq)
      enddo   !! over iq
	call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,1)
	v2r=v2r/nqdiv_c
      do ir=1,nqdiv_c
	  do j=1,n_pbmt_red
	    do i=1,n_pbmt_red
	      v4_r(i,j,ir)=v2r(ir,i,j)
	    enddo
	  enddo
      enddo   !! over ir
c ---- Dynamic part ------------------------------------------------
      do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
c --------------------------------------------------------------- 
        v2r=(0.d0,0.d0)
        do iq=1,nqdiv_c
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
          call boson_unpack_tau(tmp,n_pbmt_red,
     &                          w_tau(1,1,ind_tau,iq),npb,n_pbmt_red)
	    v2r(k11,:,:)=tmp
        enddo   !! over iq
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,1)
	  v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
	    do j=1,n_pbmt_red
	      do i=1,n_pbmt_red
	        w4_r(i,j,ir,i_tau)=v2r(ir,i,j)
	      enddo
	    enddo
        enddo   !! over ir
      enddo  !! over ind_tau
      deallocate(v2r,tmp)
	if(nproc_t/=1) then
	  ii=n_pbmt_red**2*nqdiv_c*(n_tau/2+1)
	  call dgop(w4_r,ii,'  +',comm_t)
	endif
      end