      subroutine from_tau_tau_to_tau_omega_d(n,vrt,f_tau)
c   ------- Transform F(t;t') ---> A(t;w_d) B(t;w_d) -- N elements ----
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau,2,ndim3_tau)
	complex*16, intent(out) :: vrt(n,0:n_tau,2,ndim3_omega)
      integer :: ind_tau,i,i_omega,i_tau,ind_omega,it,i1,i0,i2,nn
      real*8, allocatable :: c(:,:),cc(:,:),dd(:,:)
	complex*16, allocatable :: a_om(:,:,:,:),b_om(:,:,:,:),tmp(:,:,:),
     &	                         x(:,:),y(:,:)
      allocate(c(n,0:n_omega))
	allocate(a_om(n,0:n_omega,2,ndim3_tau))
	allocate(b_om(n,0:n_omega,2,ndim3_tau))
	allocate(cc(n,n_tau/2))
	allocate(dd(n,n_tau/2))
	allocate(x(n,0:n_omega))
	allocate(y(n,0:n_omega))
      do ind_tau=1,ndim3_tau
	  i_tau=me_t*ndim3_tau+ind_tau-1
	  do it=1,2
	    if(it==2) i_tau=n_tau-i_tau
c ------------------------ C -------------------------------------
          cc=0.d0
	    do i=1,n_bound_tau(i_tau)/2
	      i1=n_bound_tau(i_tau)-i
	      cc(:,i)=f_tau(:,i,it,ind_tau)+f_tau(:,i1,it,ind_tau)
	    enddo
	    i0=0
	    do i=n_bound_tau(i_tau)/2+1,n_tau1(i_tau)/2-1
	      i0=i0+1
	      i1=n_bound_tau(i_tau)+i0
	      i2=n_tau1(i_tau)-i0
	      cc(:,i)=f_tau(:,i1,it,ind_tau)-f_tau(:,i2,it,ind_tau)
	    enddo
c ------------------------ D -------------------------------------
          dd=0.d0
	    do i=1,n_bound_tau(i_tau)/2-1
	      i1=n_bound_tau(i_tau)-i
	      dd(:,i)=f_tau(:,i,it,ind_tau)-f_tau(:,i1,it,ind_tau)
	    enddo
	    i0=0
	    do i=n_bound_tau(i_tau)/2+1,n_tau1(i_tau)/2
	      i0=i0+1
	      i1=n_bound_tau(i_tau)+i0
	      i2=n_tau1(i_tau)-i0
	      dd(:,i)=f_tau(:,i1,it,ind_tau)+f_tau(:,i2,it,ind_tau)
	    enddo
c ------------------------------ X ------------------------------      
c ---------------------- COS contribution ----------------------
          nn=n_tau1(i_tau)/2
	    call dgemm('n','n',n,n_omega+1,nn,1.d0,cc,n,
     &               tau_omega_from_tau_tau(1,0,1,1,it,ind_tau),
     &               n_tau/2,0.d0,c,n)
          x=c
c ---------------------- SIN contribution ----------------------
	    call dgemm('n','n',n,n_omega+1,nn,1.d0,dd,n,
     &               tau_omega_from_tau_tau(1,0,2,2,it,ind_tau),
     &               n_tau/2,0.d0,c,n)
          x=x-(0.d0,1.d0)*c
c ------------------------------ Y ------------------------------      
c ---------------------- COS contribution ----------------------
	    call dgemm('n','n',n,n_omega+1,nn,1.d0,dd,n,
     &               tau_omega_from_tau_tau(1,0,1,2,it,ind_tau),
     &               n_tau/2,0.d0,c,n)
          y=c
c ---------------------- SIN contribution ----------------------
	    call dgemm('n','n',n,n_omega+1,nn,1.d0,cc,n,
     &               tau_omega_from_tau_tau(1,0,2,1,it,ind_tau),
     &               n_tau/2,0.d0,c,n)
          y=y-(0.d0,1.d0)*c
          a_om(:,:,it,ind_tau)=x+y
          b_om(:,:,it,ind_tau)=x-y
        enddo   !! over it
      enddo   !! over ind_tau
      deallocate(x,y)
      allocate(tmp(0:n_tau,2,0:n_omega))
      do i=1,n
        tmp=(0.d0,0.d0)
        do i_omega=0,n_omega
          do ind_tau=1,ndim3_tau
            i_tau=me_t*ndim3_tau+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              tmp(i_tau,1,i_omega)=a_om(i,i_omega,it,ind_tau)
              tmp(i_tau,2,i_omega)=b_om(i,i_omega,it,ind_tau)
            enddo
          enddo
        enddo
	  if(nproc_t/=1) call dgop(tmp,4*(n_tau+1)*(n_omega+1),'  +',
     & 	                           comm_t)
        do ind_omega=1,ndim3_omega
	    i_omega=me_t*ndim3_omega+ind_omega-1
	    do i_tau=0,n_tau
	      vrt(i,i_tau,:,ind_omega)=tmp(i_tau,:,i_omega)
	    enddo
	  enddo
	enddo
      deallocate(c,a_om,b_om,tmp,cc,dd)
      end