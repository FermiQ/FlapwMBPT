      subroutine vx2_stat_transf_pw(vx2_stat,vrt2_tau,ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(in) :: vx2_stat(nr_full_red,nr_full_red,
     &                                   nqdiv_c,ndim3_nu)
      complex*16, intent(out) :: vrt2_tau(nr_full_red,nr_full_red,
     &                                    nqdiv_c,2,ndim3_tau)
      integer :: it,k,ind_nu,ind_tau,k0,nk,kf0
      complex*16, allocatable :: fnu(:,:,:),ftau(:,:,:,:)
      allocate(fnu(nr_full_red,nr_full_red,ndim3_nu))
      allocate(ftau(nr_full_red,nr_full_red,2,ndim3_tau))
	vrt2_tau=(0.d0,0.d0)
	do k=1,nqdiv_c
	  k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
	  nk=n_low_bnd(kf0,ispin)
	  do ind_nu=1,ndim3_nu
	    fnu(:,:,ind_nu)=vx2_stat(:,:,k,ind_nu)
	  enddo
	  call nu_to_tau_half(fnu,ftau,nr_full_red**2)
	  do ind_tau=1,ndim3_tau
	    do it=1,2
  	      vrt2_tau(:,:,k,it,ind_tau)=ftau(:,:,it,ind_tau)
          enddo
        enddo
      enddo
      deallocate(fnu,ftau)
      end
      
      
      subroutine vx2_stat_transf_0(vx2_stat_0,vrt20)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      complex*16, intent(in) :: vx2_stat_0(maxel_red,maxel_red,natom,
     &                                     ndim3_nu)
      complex*16, intent(out) :: vrt20(maxel_red,maxel_red,natom,2,
     &                                 ndim3_tau)
      integer :: it,ind_nu,ind_tau,iatom,isort,n
      complex*16, allocatable :: fnu(:,:,:),ftau(:,:,:,:)
      vrt20=(0.d0,0.d0)
      do iatom=1,natom
        isort=is(iatom)
        n=lfunm_red(isort)
        allocate(fnu(n,n,ndim3_nu))
        allocate(ftau(n,n,2,ndim3_tau))
c --------------------------------------------------------------
	  do ind_nu=1,ndim3_nu
	    fnu(:,:,ind_nu)=vx2_stat_0(1:n,1:n,iatom,ind_nu)
	  enddo
	  call nu_to_tau_half(fnu,ftau,n*n)
	  do ind_tau=1,ndim3_tau
	    do it=1,2
  	      vrt20(1:n,1:n,iatom,it,ind_tau)=ftau(:,:,it,ind_tau)
          enddo
        enddo
        deallocate(fnu,ftau)
      enddo
      end