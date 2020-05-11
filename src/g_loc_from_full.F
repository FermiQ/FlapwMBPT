      subroutine g_loc_from_full(g_loc,g_loc_tau,z_red)
      use atom_mod
	  use heg_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      real*8, intent(out) :: g_loc_tau(nrel,maxel_red,maxel_red,2,
     &                                 ndim3_tau,nspin,natom)
      complex*16, intent(out) :: g_loc(maxel_red,maxel_red,ndim3_omega,
     &                                 nspin,natom)
      integer :: ispin,ind_k,k,ib,i,j,iatom,isort,ind0,nbnd,ind_tau,it,
     &           jb,kc,nb,jbnd,ibnd
      complex*16 :: cc
      real*8, allocatable :: gt(:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),gloc(:,:,:,:)
	  allocate(tmp1(nfun_red,nbndf))
	  allocate(tmp(nbndf_bnd,nbndf_bnd))
	  allocate(gloc(maxel_red,maxel_red,natom,2))
      allocate(gt(nbndf_bnd,nbndf_bnd,2))
c -------------------- G on TAU ------------------------------------
      g_loc_tau=0.d0
	  do ispin=1,nspin
        do ind_tau=1,ndim3_tau
          gloc=(0.d0,0.d0)
	      do ind_k=1,ndim3_k(me_k+1)
	        k=n3_mpi_k(me_k+1)+ind_k
	        nbnd=n_bnd(k,ispin)
            kc=k_c_from_a(k)
            if(kc==0) cycle
            nb=n_low_bnd(kc,ispin)
	        do jbnd=1,nb
	          jb=ind_bands_bnd(jbnd,kc,ispin)
	          do ibnd=1,nb
	            ib=ind_bands_bnd(ibnd,kc,ispin)
	            gt(ibnd,jbnd,:)=g_full(ib,jb,:,ind_tau,ind_k,ispin)
              enddo
            enddo
            do it=1,2
              call ferm_unpack_tau(tmp,gt,nb,nbndf_bnd,nbndf_bnd,it)
	          call zgemm('n','n',nfun_red,nb,nb,(1.d0,0.d0),
     &		             z_red(1,1,kc,ispin),nfun_red,tmp,nbndf_bnd,
     &                   (0.d0,0.d0),tmp1,nfun_red)
              do iatom=1,natom
	            isort=is(iatom)
	            ind0=io_lem_red(iatom)-1
                do j=1,lfunm_red(isort) 
	              do ib=1,nb
	                do i=1,lfunm_red(isort)
	                  cc=wgt_c(kc)*tmp1(ind0+i,ib)
     &                         *conjg(z_red(ind0+j,ib,kc,ispin))
	                  gloc(i,j,iatom,it)=gloc(i,j,iatom,it)+cc
                    enddo
                  enddo
                enddo 
              enddo  !! over iatom_c
            enddo  !! over it
          enddo  !! over ind_k
          do it=1,2
            call symg_red(gloc(1,1,1,it))
            do iatom=1,natom
	          isort=is(iatom)
              do j=1,lfunm_red(isort) 
	            do i=1,lfunm_red(isort)
	              g_loc_tau(1,i,j,it,ind_tau,ispin,iatom)=
     &	                                       gloc(i,j,iatom,it)
	            enddo
              enddo
              if(irel==2) then
                do j=1,lfunm_red(isort)
                  do i=1,lfunm_red(isort)
                    g_loc_tau(2,i,j,it,ind_tau,ispin,iatom)=
     &                         imag(gloc(i,j,iatom,it))
                  enddo
                enddo
              endif
            enddo
          enddo   !! over it
	    enddo   !! over ind_tau
	  enddo   !! over ispin
	  if(nproc_k/=1) call dgop(g_loc_tau,
     &	                     2*nrel*maxel_red**2*ndim3_tau*nspin*natom,
     &                         '  +',comm_k)
c -------------------- G on OMEGA ------------------------------------
      do iatom=1,natom
	    isort=is(iatom)
        if(lfunm_red(isort)==0) cycle
	    do ispin=1,nspin
	      if(irel/=2) then
            call tau_to_omega_abb(g_loc(1,1,1,ispin,iatom),maxel_red,
     &                            g_loc_tau(1,1,1,1,1,ispin,iatom),
     &                            maxel_red,lfunm_red(isort))
          else if(irel==2) then
            call tau_to_omega_abbc(g_loc(1,1,1,ispin,iatom),maxel_red,
     &                             g_loc_tau(1,1,1,1,1,ispin,iatom),
     &                             maxel_red,lfunm_red(isort))
          endif
        enddo
      enddo
	  deallocate(tmp,tmp1,gloc)
      end
