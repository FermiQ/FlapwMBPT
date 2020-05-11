      subroutine psi_psi_pw_red
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer :: k,ispin,k0,ibnd,ib,ind_k,kf,ka,ka0
	  real*8, allocatable :: fifim(:,:,:,:,:)
	  complex*16, allocatable :: zw(:,:,:,:),aw(:,:,:,:)
c ------------------------------------------------------------------
      allocate(fifim(maxel,maxel,n_pbmtm,nsort,nspin))
c ---------- Interaction matrices -----------------------------------
      do ispin=1,nspin
        call fi0_full(fifim(1,1,1,1,ispin),ispin,ispin)
      enddo
c ----------------------------------------------------------
	  allocate(zw(nfun,nbndf_bnd,nqdiv_c,nspin))
	  allocate(aw(nbasmpw,nbndf_bnd,nqdiv_c,nspin))
      aw=(0.d0,0.d0)
      zw=(0.d0,0.d0)
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me_k+1)
	      kf=n3_mpi_k(me_k+1)+ind_k
          k=k_c_from_a(kf)
          if(k==0) cycle
	      do ibnd=1,n_low_bnd(k,ispin)
	        ib=ind_bands_bnd(ibnd,k,ispin)
	        zw(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
	        aw(:,ibnd,k,ispin)=ev_bnd(:,ib,ind_k,ispin)
	      enddo
	    enddo
        do k=npnt_c+1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(zw(1,1,k,ispin),zw(1,1,k0,ispin),
     &               k_group_c(k),pnt_c(1,k),nfun,
     &               n_low_bnd(ka0,ispin),maxb,io_lem,lmb,indbasa,limlb,
     &               0)
	      call sym_a(aw(1,1,k,ispin),aw(1,1,k0,ispin),ka,
     &               k_group_c(k),pnt_c(1,k),nbasmpw,
     &               n_low_bnd(ka0,ispin),nbask(ka0),
     &               indgb(1,ka),iplf_bk(1,ka0))
        enddo
      enddo
	  if(nproc_k/=1) then
	    call dgop(zw,2*nfun*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
	    call dgop(aw,2*nbasmpw*nbndf_bnd*nqdiv_c*nspin,'  +',comm_k)
	  endif
	  if(vrt_x_appr(2:3)/='00') then
	    if(allocated(psi_psi_pw)) deallocate(psi_psi_pw)
	    allocate(psi_psi_pw(nbndf_bnd,nbndf_bnd,nqdiv_c,nc_line,nspin,
     &                      nspin))
	    call psi_psi_pw_prepare(zw,aw)
	  endif
	  deallocate(fifim,aw,zw)
      end
