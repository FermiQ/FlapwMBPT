      subroutine sig_omega0_solid(n,sigc,ind_k,ispin,sigc0,zren)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer, intent(in) :: n,ind_k,ispin
	  real*8, intent(in) :: sigc(n,n,2,ndim3_omega)
	  complex*16, intent(out) :: sigc0(n,n),zren(n,n)
	  integer :: i,j,ind_omega,i_omega
	  complex*16, allocatable :: tmp(:,:,:),tmp0(:,:)
	  allocate(tmp(n,n,2))
	  allocate(tmp0(n,0:n_omega))
	  tmp=(0.d0,0.d0)
	  tmp0=(0.d0,0.d0)
	  do ind_omega=1,ndim3_omega
	    i_omega=me_t*ndim3_omega+ind_omega-1
	    do j=1,n
	      tmp0(j,i_omega)=dcmplx(sigc(j,j,1,ind_omega),
     &                           sigc(j,j,2,ind_omega))
        enddo
	    if(i_omega<=1) then
	      i=i_omega+1
          call ferm_unpack_omega(tmp(1,1,i),sigc(1,1,1,ind_omega),n,n,n)
        endif
      enddo   !! over ind_omega
	  if(nproc_t/=1) then
	    call dgop(tmp0,2*n*(n_omega+1),'  +',comm_t)
	    call dgop(tmp,4*n*n,'  +',comm_t)
	  endif
	  call sigma_0_z_ren(n,n,tmp,tmp0)
	  sigc0=tmp(:,:,1)
	  zren=tmp(:,:,2)
	  deallocate(tmp,tmp0)
      end
