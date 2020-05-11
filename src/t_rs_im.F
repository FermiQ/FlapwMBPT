      subroutine t_rs_im(it,ll,t_pw,t_rs,phf)
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
     &                               2),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(nr_full_red,n_pbmt_red,nqdiv_c,ll)
	  integer :: ir,j,i,k11,il,iq1,iq10,np,ka,ka0,i0,ind
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:)
      allocate(v2r(nqdiv_c,nr_full_red,n_pbmt_red))
      allocate(tmp(nr_full_red,n_pbmt_red))
      do il=1,ll
        v2r=(0.d0,0.d0)
	    do iq1=1,nqdiv_c
	      call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,k11)
          iq10=i_kref_c(iq1)
          ka0=k_a_from_c(iq10)
          ka=k_a_from_c(iq1)
	      np=nplwgw_red(ka0)
          tmp=(0.d0,0.d0)
          do i=1,np
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            tmp(ind,:)=t_pw(n_pbmt_red+i,1:n_pbmt_red,iq1,il,3-it)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              tmp,1)
          do i=1,nr_full_red
            v2r(k11,i,:)=phf(i,iq1)*tmp(i,:)
          enddo
        enddo
        v2r=v2r/sqrt(amega)
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &            v2r,-1)
	    v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
          do i=1,nr_full_red
            do j=1,n_pbmt_red
              t_rs(i,j,ir,il)=v2r(ir,i,j)
            enddo
          enddo
        enddo   !! over ir
      enddo  !! over il
      deallocate(v2r,tmp)
      end
