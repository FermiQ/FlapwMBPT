      subroutine from_tau_tau_to_tau_nu_d(n,vrt,f_tau)
c   ------- Transform F(t;t') ---> A(t;v) B(t;v) -- N elements ----
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau,2,ndim3_tau)
	complex*16, intent(out) :: vrt(n,0:n_tau,2,ndim3_nu)
	complex*16, allocatable :: a_nu(:,:,:,:),b_nu(:,:,:,:),tmp(:,:,:)
      integer :: ind_tau,i,i_nu,i_tau,ind_nu,it
      real*8, allocatable :: c(:,:)
      allocate(c(n,0:n_nu))
	allocate(a_nu(n,0:n_nu,2,ndim3_tau))
	allocate(b_nu(n,0:n_nu,2,ndim3_tau))
      do ind_tau=1,ndim3_tau
        do it=1,2
c ------------------------------ A ------------------------------      
c ---------------------- COS contribution ----------------------
	    call dgemm('n','n',n,n_nu+1,n_tau+1,1.d0,
     &	             f_tau(1,0,it,ind_tau),n,
     &               tau_nu_from_tau_tau(0,0,1,1,it,ind_tau),n_tau+1,
     &               0.d0,c,n)
          a_nu(:,:,it,ind_tau)=c      
c ---------------------- SIN contribution ----------------------
	    call dgemm('n','n',n,n_nu+1,n_tau+1,1.d0,
     &	             f_tau(1,0,it,ind_tau),n,
     &               tau_nu_from_tau_tau(0,0,2,1,it,ind_tau),n_tau+1,
     &               0.d0,c,n)
          a_nu(:,:,it,ind_tau)=a_nu(:,:,it,ind_tau)+(0.d0,1.d0)*c
c ------------------------------ B ------------------------------      
c ---------------------- COS contribution ----------------------
	    call dgemm('n','n',n,n_nu+1,n_tau+1,1.d0,
     &	             f_tau(1,0,it,ind_tau),n,
     &               tau_nu_from_tau_tau(0,0,1,2,it,ind_tau),n_tau+1,
     &               0.d0,c,n)
          b_nu(:,:,it,ind_tau)=c      
c ---------------------- SIN contribution ----------------------
	    call dgemm('n','n',n,n_nu+1,n_tau+1,1.d0,
     &	             f_tau(1,0,it,ind_tau),n,
     &               tau_nu_from_tau_tau(0,0,2,2,it,ind_tau),n_tau+1,
     &               0.d0,c,n)
          b_nu(:,:,it,ind_tau)=b_nu(:,:,it,ind_tau)-(0.d0,1.d0)*c
        enddo   !! over it
      enddo   !! over ind_tau
      allocate(tmp(0:n_tau,2,0:n_nu))
      do i=1,n
        tmp=(0.d0,0.d0)
        do i_nu=0,n_nu
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              tmp(i_tau,1,i_nu)=a_nu(i,i_nu,it,ind_tau)
              tmp(i_tau,2,i_nu)=b_nu(i,i_nu,it,ind_tau)
            enddo
          enddo
        enddo
	  if(nproc_t/=1) call dgop(tmp,2*(n_tau+1)*(n_nu+1),'  +',
     & 	                           comm_t)
        do ind_nu=1,ndim3_nu
	    i_nu=me_t*ndim3_nu+ind_nu-1
	    do i_tau=0,n_tau
	      vrt(i,i_tau,:,ind_nu)=tmp(i_tau,:,i_nu)
	    enddo
	  enddo
	enddo
      deallocate(c,a_nu,b_nu,tmp)
      end