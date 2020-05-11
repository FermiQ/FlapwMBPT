      subroutine v4_for_bnd_ii(vii,v2_q,phf)
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
	  real*8, intent(out) :: vii(nr_red_red,nr_full_red,nqdiv_c)
	  integer :: ir,iq0,j,i,mx2,iq,k11,ka0,i0,ka,np,ind
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:),t1(:,:)
c ----------------------------------------------------------------
      vii=0.d0
      mx2=maxel_red**2
      allocate(v2r(nqdiv_c,nr_red_red,nr_full_red))
      allocate(t1(nr_full_red,nr_full_red))
      allocate(tmp(nr_full_red,nplw_gw_red))
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
          tmp(ind,1:np)=v2_q(n_pbmt_red+1:n_pbmt_red+np,n_pbmt_red+i,iq)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),np,tmp,-1)
        do i=1,nr_full_red
          tmp(i,:)=conjg(phf(i,iq))*tmp(i,:)
        enddo
        t1=(0.d0,0.d0)
        do i=1,np
          i0=indpw_gw_red(i,ka)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t1(ind,:)=tmp(:,i)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,t1,
     &            1)
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
        do i0=1,nr_red_red
          i=r_pnt_ind_red(i0)
          v2r(k11,i0,:)=phf(i,iq)*t1(i,:)
        enddo
      enddo   !! over iq
      v2r=v2r/amega
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nr_full_red*nr_red_red,
     &          v2r,1)
	  v2r=v2r/nqdiv_c
      do ir=1,nqdiv_c
	    do j=1,nr_full_red
	      do i=1,nr_red_red
            vii(i,j,ir)=v2r(ir,i,j)
          enddo  !! over i
        enddo   !! over j
      enddo  !! over ir
      deallocate(v2r,tmp,t1)
      end
