      subroutine vx2_stat_transf(vx2_stat,vrt2_tau,vrt2_nu,z_red,a_red,
     &                           phf,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: vx2_stat(nrs_red,nrs_red,nqdiv_c,
     &                                   ndim3_nu),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      complex*16, intent(out) :: vrt2_tau(nbndf_bnd,nrs_red,nqdiv_c,2,
     &                                    ndim3_tau),
     &                           vrt2_nu(nbndf_bnd,nrs_red,nqdiv_c,
     &                                   ndim3_nu)
      integer :: it,k,ind_nu,ind_tau,k0,nk,kf0
      complex*16, allocatable :: tm(:,:),fnu(:,:,:),ftau(:,:,:,:)
      allocate(tm(nrs_red,nbndf_bnd))
      allocate(fnu(nbndf_bnd,nrs_red,ndim3_nu))
      allocate(ftau(nbndf_bnd,nrs_red,2,ndim3_tau))
	  vrt2_tau=(0.d0,0.d0)
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	    nk=n_low_bnd(kf0,ispin)
c --------------------------------------------------------------
	    do ind_nu=1,ndim3_nu
c --- Transform VRT_(tL)r_(t'L')r' to VRT_lambda_(tL')r' -----------
            call vrt_r_r1_to_bnd_r1(k,nk,z_red(1,1,k),a_red(1,1,k),phf,
     &                              vx2_stat(1,1,k,ind_nu),
     &                              vrt2_nu(1,1,k,ind_nu))
	      fnu(:,:,ind_nu)=vrt2_nu(:,:,k,ind_nu)
	    enddo
	    call nu_to_tau_half(fnu,ftau,nbndf_bnd*nrs_red)
	    do ind_tau=1,ndim3_tau
	      do it=1,2
  	        vrt2_tau(:,:,k,it,ind_tau)=ftau(:,:,it,ind_tau)
          enddo
        enddo
      enddo
      deallocate(tm,fnu,ftau)
      end
