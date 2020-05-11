      subroutine t_rs_mm(it,ll,t_pw,t_rs)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: it,ll
      complex*16, intent(in) :: t_pw(n_pbtot_red,n_pbtot_red,nqdiv_c,ll,
     &                               2)
      complex*16, intent(out) :: t_rs(n_pbmt_red,n_pbmt_red,nqdiv_c,ll)
	integer :: ir,il,ind,j,i,iq1
	complex*16, allocatable :: v2r(:,:,:)
c ----------------------------------------------------------------
      t_rs=(0.d0,0.d0)
      allocate(v2r(nqdiv_c,n_pbmt_red,n_pbmt_red))
      do il=1,ll
c --------------------------------------------------------------- 
        v2r=(0.d0,0.d0)
	  do iq1=1,nqdiv_c
          call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,ind)
	    v2r(ind,:,:)=t_pw(1:n_pbmt_red,1:n_pbmt_red,iq1,il,3-it)
        enddo   !! over iq1
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,-1)
	  v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
          do j=1,n_pbmt_red
            do i=1,n_pbmt_red
              t_rs(i,j,ir,il)=v2r(ir,i,j)
            enddo
          enddo
        enddo   !! over ind_r
      enddo  !! over in
      deallocate(v2r)
      end