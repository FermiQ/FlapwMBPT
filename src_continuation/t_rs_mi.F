      subroutine t_rs_mi(it,ll,t_pw,t_rs,iq,phf)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: it,iq,ll
      complex*16, intent(in) :: t_pw(n_pbtot_red,n_pbtot_red,nqdiv_c,ll,
     &                               2),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(n_pbmt_red,nr_full_red,nqdiv_c,ll)
	  integer :: ir,j,i,k11,il,iq1,iq10,iq2,iq20,ka,ka0,np,i0,ind
      real*8 :: v(3)
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:)
c ----------------------------------------------------------------
      allocate(v2r(nqdiv_c,n_pbmt_red,nr_full_red))
      allocate(tmp(nr_full_red,n_pbmt_red))
      do il=1,ll
        v2r=(0.d0,0.d0)
	    do iq1=1,nqdiv_c
          iq10=i_kref_c(iq1)
          v=pnt_c(:,iq1)-pnt_c(:,iq)
          call zone1_number(v,rb0_c,ndiv_c,iq2)  ! iq2 - index in BZ_1
          iq2=index_k1_c(iq2)  ! Now iq2 - index in BZ_0
	      iq20=i_kref_c(iq2)
          ka0=k_a_from_c(iq20)
          ka=k_a_from_c(iq2)
	      np=nplwgw_red(ka0)
          tmp=(0.d0,0.d0)
          do i=1,np
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            tmp(ind,:)=t_pw(1:n_pbmt_red,n_pbmt_red+i,iq1,il,3-it)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              tmp,-1)
	      call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,k11)
          do i=1,nr_full_red
            v2r(k11,:,i)=conjg(phf(i,iq2))*tmp(i,:)
          enddo
        enddo   !! over iq1
        v2r=v2r/sqrt(amega)
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &            v2r,-1)
	    v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            do i=1,n_pbmt_red
              t_rs(i,j,ir,il)=v2r(ir,i,j)
            enddo
          enddo
        enddo   !! over ind_r
      enddo  !! over in
      deallocate(v2r,tmp)
      end
