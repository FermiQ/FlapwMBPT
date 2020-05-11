      subroutine t_rs_ii(it,ll,t_pw,t_rs,phf,iq)
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
	  complex*16, intent(out) :: t_rs(nr_full_red,nr_full_red,nqdiv_c,
     &                                ll)
      integer :: ir,np,k11,il,iq1,iq10,iq2,j,i,iq20,ka,ka0,kqa0,kqa,np1,
     &           i0,ind
      real*8 :: v(3)
	  complex*16, allocatable :: v2r(:,:,:),t(:,:),t1(:,:)
c ----------------------------------------------------------------
      allocate(v2r(nqdiv_c,nr_full_red,nr_full_red))
      allocate(t(nr_full_red,nplw_gw_red))
      allocate(t1(nr_full_red,nr_full_red))
      do il=1,ll
        v2r=(0.d0,0.d0)
	    do iq1=1,nqdiv_c
          iq10=i_kref_c(iq1)
          ka0=k_a_from_c(iq10)
          ka=k_a_from_c(iq1)
	      np=nplwgw_red(ka0)
          v=pnt_c(:,iq1)-pnt_c(:,iq)
          call zone1_number(v,rb0_c,ndiv_c,iq2)  ! iq2 - index in BZ_1
          iq2=index_k1_c(iq2)  ! Now iq2 - index in BZ_0
	      iq20=i_kref_c(iq2)
          kqa0=k_a_from_c(iq20)
          kqa=k_a_from_c(iq2)
	      np1=nplwgw_red(kqa0)
	      call zone1_number(pnt_c(1,iq1),rb0_c,ndiv_c,k11)
          t=(0.d0,0.d0)
          do i=1,np1
            i0=indpw_gw_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t(ind,1:np)=t_pw(n_pbmt_red+1:n_pbmt_red+np,n_pbmt_red+i,
     &                       iq1,il,3-it)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),np,t,-1)
          do i=1,nr_full_red
            t(i,1:np)=conjg(phf(i,iq2))*t(i,1:np)
          enddo
          t1=(0.d0,0.d0)
          do i=1,np
            i0=indpw_gw_red(i,ka)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(ind,:)=t(:,i)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,
     &              t1,1)
          do i=1,nr_full_red
            v2r(k11,i,:)=phf(i,iq1)*t1(i,:)
          enddo
        enddo   !! over iq1
        v2r=v2r/amega
	    call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nr_full_red**2,v2r,-1)
	    v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            do i=1,nr_full_red
              t_rs(i,j,ir,il)=v2r(ir,i,j)
            enddo
          enddo
        enddo   !! over ir
      enddo  !! over il
      deallocate(v2r,t,t1)
      end
