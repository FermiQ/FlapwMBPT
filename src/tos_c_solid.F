      subroutine tos_c_solid(x,ttos)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
	  real*8, intent(in) :: x
	  real*8, intent(inout) :: ttos
	  integer :: k,ispin,ndim,i_omega1,i_omega,i,ind_k
	  real*8 :: const,st,w_n,ttos_c
	  real*8, allocatable :: gc_diag(:,:),sigc(:,:,:,:)
	  complex*16, allocatable :: gc_omega(:,:),gx_omega(:)
	  const=2.d0/dfloat(nspin*nrel)
	  ttos_c=0.d0
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
	    do ispin=1,nspin
	      ndim=n_bnd(k,ispin)
          allocate(sigc(ndim,ndim,2,ndim3_omega))
          call tau_to_omega_baa(sigc,ndim,
     &                          sig_c_tau(1,1,1,1,ind_k,ispin),nbndf,
     &                          ndim)
	      allocate(gc_omega(ndim,ndim),gx_omega(ndim),
     &             gc_diag(ndim,0:n_omega))
          gc_diag=0.d0
	      do i_omega1=1,ndim3_omega
	        i_omega=me_t*ndim3_omega+i_omega1-1
	        w_n=w_omega(i_omega)
    	    call g_x_omega(ispin,k,gx_omega,w_n,ndim,x)
c ------ We temporarily place SIGMA_C into GC_OMEGA ----------------
            call ferm_unpack_omega(gc_omega,sigc(1,1,1,i_omega1),ndim,
     &                             ndim,ndim)
            call g_c_omega(0,gx_omega,gc_omega,ndim,ndim)
c ----------------------------------------------------------------------
	        do i=1,ndim
	          gc_diag(i,i_omega)=real(gc_omega(i,i))
	        enddo
	      enddo  !! over i_omega1
	      deallocate(sigc)
	      if(nproc_t/=1) call dgop(gc_diag,ndim*(n_omega+1),'  +',
     &	                             comm_t)
	      st=0.d0
	      do i=1,ndim
            do i_omega=0,n_omega
              st=st+tau_from_omega(i_omega,0,1)*gc_diag(i,i_omega)
            enddo
          enddo  !! over idim
	      st=-st/2.d0
	      deallocate(gc_diag)
	      deallocate(gx_omega)
	      deallocate(gc_omega)
          ttos_c=ttos_c-wgt(k)*st*const
	    enddo  !! over ispin
	  enddo  !! over ind_k
      if(nproc_k/=1) call DGOP(ttos_c,1,'  +',comm_k)
      if(maswrk) then
        write(iun,'(1x,f13.8,3(1x,f13.8))')x,ttos,ttos_c,ttos+ttos_c
      endif
	  ttos=ttos+ttos_c
      end
      
      
      subroutine tos_qp_solid(x,ttos)
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      use vertex_mod
      implicit none
	  real*8, intent(in) :: x
	  real*8, intent(inout) :: ttos
	  integer :: me1,ispin,ind_k,k,i,n
	  real*8 :: const,st,de,ffd,fermi_dirac
	  me1=me+1
	  const=2.d0/dfloat(nspin*nrel)
	  call bands_qp(x)
      do ispin=1,nspin
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n=n_bnd(k,ispin)
          st=0.d0
	      do i=1,n
            de=e_qp(i,k,ispin)-x
            ffd=fermi_dirac(de)
            st=st+ffd*qp_spectr_weight(i,ind_k,ispin)
          enddo
          ttos=ttos+wgt(k)*st*const
        enddo  !!! over k
      enddo  !!! over ispin
      if(nproc_k/=1) call DGOP(ttos,1,'  +',comm_k)
      if(maswrk) write(iun,'(3(1x,f13.8))')x,ttos
      end
