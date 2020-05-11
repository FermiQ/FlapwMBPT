      subroutine aux_full(fif_full,z_full,a_full)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(out) :: fif_full(nrel,maxel**2,n_pbmtm_red,nsort,
     &                                nspin)
      complex*16, intent(out) :: z_full(nfun,nbndf_bnd,nqdiv_c,nspin),
     &                           a_full(nbasmpw,nbndf_bnd,nqdiv_c,
     &                                  nspin)
      integer :: k,ispin,i,i1,isort,ij,k0,ka,ka0,ind_k,ibnd,ib,irl
      real*8 :: pi2
      real*8, allocatable :: fifim(:,:,:,:,:)
      pi2=pi+pi
      fif_full=0.d0
c ---------------------------------------------------------------            
      allocate(fifim(nrel,maxel,maxel,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red_1(fifim,ispin,ispin)
        do isort=1,nsort
          ij=0
          do i1=1,lfunm(isort)
            do i=1,lfunm(isort)
              ij=ij+1
              do irl=1,nrel
                fif_full(irl,ij,:,isort,ispin)=fifim(irl,i,i1,:,isort)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
c -------- Z_FULL, A_FULL for K=1,nqdiv_c -------------------------------
      z_full=(0.d0,0.d0)
      a_full=(0.d0,0.d0)
	  do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          ka=n3_mpi_k(me_k+1)+ind_k
          k=k_c_from_a(ka)
          if(k==0) cycle
          do ibnd=1,n_low_bnd(ka,ispin)
            ib=ind_bands_bnd(ibnd,ka,ispin)
	        z_full(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
	        a_full(:,ibnd,k,ispin)=ev_bnd(:,ib,ind_k,ispin)
          enddo
	    enddo
        if(nproc_k/=1) then
          call dgop(z_full(1,1,1,ispin),2*nfun*nbndf_bnd*npnt_c,'  +',
     &              comm_k)
          call dgop(a_full(1,1,1,ispin),2*nbasmpw*nbndf_bnd*npnt_c,
     &              '  +',comm_k)
        endif
        do k=npnt_c+1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(z_full(1,1,k,ispin),z_full(1,1,k0,ispin),
     &               k_group_c(k),pnt_c(1,k),nfun,
     &               n_low_bnd(ka0,ispin),maxb,io_lem,lmb,indbasa,limlb,
     &               0)
	      call sym_a(a_full(1,1,k,ispin),a_full(1,1,k0,ispin),ka,
     &               k_group_c(k),pnt_c(1,k),nbasmpw,
     &               n_low_bnd(ka0,ispin),nbask(ka0),indgb(1,ka),
     &               iplf_bk(1,ka0))
        enddo
      enddo
      call timel('**** AUX_FULL finished *************')
      end
      
      
      subroutine aux_full_loc(fif_full,z_full)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(out) :: fif_full(maxel**2,n_pbmtm_red,nsort,nspin)
      complex*16, intent(out) :: z_full(nfun,nbndf_bnd,nqdiv_c,nspin)
      integer :: k,ispin,i,i1,isort,ij,k0,ka,ka0,ind_k,ibnd,ib
      real*8 :: pi2
      real*8, allocatable :: fifim(:,:,:,:)
      pi2=pi+pi
      fif_full=0.d0
c ---------------------------------------------------------------            
      allocate(fifim(maxel,maxel,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red_1(fifim,ispin,ispin)
        do isort=1,nsort
          ij=0
          do i1=1,lfunm(isort)
            do i=1,lfunm(isort)
              ij=ij+1
              fif_full(ij,:,isort,ispin)=fifim(i,i1,:,isort)
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
c -------- Z_FULL, for K=1,nqdiv_c -------------------------------
      z_full=(0.d0,0.d0)
	do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          ka=n3_mpi_k(me_k+1)+ind_k
          k=k_c_from_a(ka)
          if(k==0) cycle
          do ibnd=1,n_low_bnd(ka,ispin)
            ib=ind_bands_bnd(ibnd,ka,ispin)
	      z_full(:,ibnd,k,ispin)=z_bnd(:,ib,ind_k,ispin)
          enddo
	  enddo
        if(nproc_k/=1) call dgop(z_full(1,1,1,ispin),
     &                           2*nfun*nbndf_bnd*npnt_c,'  +',
     &                           comm_k)
        do k=npnt_c+1,nqdiv_c
	    k0=i_kref_c(k)
          ka=k_a_from_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(z_full(1,1,k,ispin),z_full(1,1,k0,ispin),
     &               k_group_c(k),pnt_c(1,k),nfun,
     &               n_low_bnd(ka0,ispin),maxb,io_lem,lmb,indbasa,limlb,
     &               0)
        enddo
      enddo
      call timel('**** AUX_FULL_LOC finished *********')
      end
