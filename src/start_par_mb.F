      subroutine start_par_mb
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      implicit none
#ifdef MPI
      include 'mpif.h'
	  integer :: ierr,color
#endif
c ---- Memory distribution for Double Parallelization  -----------------------
	  if(nproc/=nproc_t*nproc_k*nproc_b) then
	    if(me==0) write(iun,*)' NPROC /= nproc_t*NPROC_K*nproc_b'
	    call ending
	  endif
	  if(mod(n_tau/2+1,nproc_t)/=0) then
	    if(me==0) write(iun,*)' N_TAU does not match nproc_t'
	    call ending
	  endif
	  if(mod(n_omega+1,nproc_t)/=0) then
	    if(me==0) write(iun,*)' N_OMEGA does not match nproc_t'
	    call ending
	  endif
	  if(mod(n_nu+1,nproc_t)/=0) then
	    if(me==0) write(iun,*)' N_NU does not match nproc_t'
	    call ending
	  endif
	  ndim3_tau=(n_tau/2+1)/nproc_t
	  ndim3_omega=(n_omega+1)/nproc_t
	  ndim3_nu=(n_nu+1)/nproc_t
	  ndim3_tn=max(ndim3_tau,ndim3_nu)
	  allocate(ndim_tau(nproc))
	  allocate(ndim_istar(nproc))
	  allocate(n_mpi_tau(nproc))
	  allocate(n_mpi_omega(nproc))
	  allocate(n_mpi_nu(nproc))
	  allocate(n_mpi_istar(nproc))
c ---- Memory distribution for Tau-mesh ----------------------------
	  call size_shift_par(n_tau/2+1,nproc,ndim_tau,n_mpi_tau)
c ------------ Communicators ---------------------------------------
      if(mpi_pref=='tkb') then
	    me_t=mod(me,nproc_t)
	    me_k=mod(me/nproc_t,nproc_k)
	    me_b=mod(me/(nproc_t*nproc_k),nproc_b)
      else if(mpi_pref=='btk') then
        me_b=mod(me,nproc_b)
        me_t=mod(me/nproc_b,nproc_t)
        me_k=mod(me/(nproc_b*nproc_t),nproc_k)
      else if(mpi_pref=='tbk') then
        me_t=mod(me,nproc_t)
        me_b=mod(me/nproc_t,nproc_b)
        me_k=mod(me/(nproc_b*nproc_t),nproc_k)
      else if(mpi_pref=='ktb') then
        me_k=mod(me,nproc_k)
        me_t=mod(me/nproc_k,nproc_t)
        me_b=mod(me/(nproc_k*nproc_t),nproc_b)
      else if(mpi_pref=='bkt') then
        me_b=mod(me,nproc_b)
        me_k=mod(me/nproc_b,nproc_k)
        me_t=mod(me/(nproc_b*nproc_k),nproc_t)
      else if(mpi_pref=='kbt') then
        me_k=mod(me,nproc_k)
        me_b=mod(me/nproc_k,nproc_b)
        me_t=mod(me/(nproc_b*nproc_k),nproc_t)
	  endif
      me_t_b=nproc_t*me_b+me_t
      me_t_k=nproc_t*me_k+me_t
      me_k_b=nproc_k*me_b+me_k
#ifdef MPI
      call output_map_mpi
! -- Processes with the same me_k
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,me_k,me,comm_t_b,ierr)
! -- Processes with the same me_t
      call MPI_COMM_SPLIT(MPI_COMM_WORLD,me_t,me,comm_k_b,ierr)
! -- Processes with the same me_b + me_k
      color=1000000*(me_b+1)+1000+me_k
	  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,me,comm_t,ierr)
! -- Processes with the same me_b + me_t
      color=1000000*(me_b+1)+1000+me_t
	  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,me,comm_k,ierr)
! -- Processes with the same me_k + me_t
      color=1000000*(me_k+1)+1000+me_t
	  call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,me,comm_b,ierr)
#endif	
	  allocate(ndim_omega(nproc_k))
	  allocate(ndim_nu(nproc_k))
c ---- Memory distribution for Omega-mesh ----------------------------
	  call size_shift_par(n_omega+1,nproc_k,ndim_omega,n_mpi_omega)
c ---- Memory distribution for Nu-mesh -------------------------------
	  call size_shift_par(n_nu+1,nproc_k,ndim_nu,n_mpi_nu)
	  call map_pbr_scalapack
      end
