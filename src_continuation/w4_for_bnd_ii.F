      subroutine w4_for_bnd_ii(vii,wii,w_tau,v2_q,phf,nrr3,indrr3)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: nrr3
	  real*8, intent(in) :: w_tau(n_pbtot_red,n_pbtot_red,
     &                            ndim3_tau,nqdiv_c)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      integer, intent(out) :: indrr3(nr_full_red,nr_full_red,nqdiv_c)
	  real*8, intent(out) :: vii(nrr3),wii(nrr3,0:n_tau/2)
	  integer :: ind_tau,ir,iq0,j,i,i_tau,mx2,iq,k11,ii,ka0,i0,ka,np,
     &           ind,ig,ir0,j0
      real*8 :: v(3),v1(3)
	  complex*16, allocatable :: v2r(:,:,:),tmp(:,:),tt(:,:),t1(:,:)
c ----------------------------------------------------------------
      wii=0.d0
      vii=0.d0
      mx2=maxel_red**2
      allocate(v2r(nqdiv_c,nr_red_red,nr_full_red))
      allocate(t1(nr_full_red,nr_full_red))
      allocate(tt(n_pbtot_red,n_pbtot_red))
      allocate(tmp(nr_full_red,nplw_gw_red))
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
	  ii=0
      do ir=1,nqdiv_c
	    do j=1,nr_full_red
	      do i=1,nr_red_red
	        ii=ii+1
            vii(ii)=v2r(ir,i,j)
          enddo  !! over i
        enddo   !! over j
      enddo  !! over ir
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
            tmp(ind,1:np)=tt(n_pbmt_red+1:n_pbmt_red+np,n_pbmt_red+i)
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
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,
     &              t1,1)
	      call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
          do i0=1,nr_red_red
            i=r_pnt_ind_red(i0)
            v2r(k11,i0,:)=phf(i,iq)*t1(i,:)
          enddo
        enddo   !! over iq
        v2r=v2r/amega
        call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nr_full_red*nr_red_red,
     &            v2r,1)
	    v2r=v2r/nqdiv_c
c ------------------------------------------------------------------
        ii=0
        do ir=1,nqdiv_c
	      do j=1,nr_full_red
	        do i=1,nr_red_red
	          ii=ii+1
              wii(ii,i_tau)=v2r(ir,i,j)
            enddo  !! over i
          enddo   !! over j
        enddo  !! over ind_r
      enddo  !! over ind_tau
      deallocate(tt,v2r,tmp,t1)
	  if(nproc_t/=1) then
	    ii=nrr3*(n_tau/2+1)
	    call dgop(wii,ii,'  +',comm_t)
	  endif
c ------------ Indexes for optimization ----------------------------
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
            call zone1_number(v1,gbas,ndiv_c,ir0)
            j0=ip_r_red(j,ig)
            indrr3(i,j,ir)=nr_red_red*nr_full_red*(ir0-1)
     &                    +nr_red_red*(j0-1)+i0
          enddo
        enddo
      enddo
      end
