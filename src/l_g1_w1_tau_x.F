      subroutine l_g1_w1_tau_x(eloc,xloc,n,np,lx,i_nu,ffb)
      use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: n,np,i_nu
      real*8, intent(in) :: eloc(n),ffb(n,n,np)
      complex*16, intent(in) :: xloc(n,n)
      complex*16, intent(out) :: lx(n,n,np,0:n_tau)
      integer :: i,j,ii,i_tau,l,k
      real*8 :: de,gt,green0,gt0,de1
      complex*16, allocatable :: ll(:,:),ff(:,:),xx(:,:),y(:,:,:)
      allocate(xx(n,n))
      allocate(ll(n,n))
      allocate(y(n,n,np))
      allocate(ff(n,n))
      do ii=1,np
        ff=ffb(:,:,ii)
	  call zgemm('c','n',n,n,n,(1.d0,0.d0),xloc,n,ff,n,(0.d0,0.d0),xx,
     &	           n)
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),xx,n,xloc,n,(0.d0,0.d0),
     &	           y(1,1,ii),n)
      enddo
      lx=(0.d0,0.d0)
      do i_tau=0,n_tau
        do j=1,n
          de=eloc(j)-chem_pot
	    gt=green0(de,tau_mesh(i_tau))
	    do i=1,n
	      de1=eloc(j)-eloc(i)
	      if(i_nu==0.and.abs(de1)<1.d-10) then
	        gt0=green0(de,betta_t)
	        ll(i,j)=gt*(betta_t*gt0+tau_mesh(i_tau))/2.d0
	      else
	        ll(i,j)=-gt/dcmplx(de1,w_nu(i_nu))
	      endif
	    enddo
	  enddo
	  do ii=1,np
	    do j=1,n
	      do i=1,n
	        do l=1,n
	          do k=1,n
	            lx(i,j,ii,i_tau)=lx(i,j,ii,i_tau)+xloc(i,k)*y(k,l,ii)
     &	                                           *ll(k,l)
     &                                             *conjg(xloc(j,l))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
	deallocate(ll,xx,y,ff)
      end
      
      
      
      subroutine l1_g1_w1_tau_x(eloc,xloc,n,np,lx,i_nu,lam)
      use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: n,np,i_nu
      real*8, intent(in) :: eloc(n)
      complex*16, intent(in) :: lam(n,n,np),xloc(n,n)
      complex*16, intent(out) :: lx(n,n,np,0:n_tau)
      integer :: i,j,ii,i_tau,l,k
      real*8 :: de,gt,green0,gt0,de1
      complex*16, allocatable :: ll(:,:),xx(:,:),y(:,:,:)
      allocate(xx(n,n))
      allocate(ll(n,n))
      allocate(y(n,n,np))
      y=(0.d0,0.d0)
      do ii=1,np
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),lam(1,1,ii),n,xloc,n,
     &	           (0.d0,0.d0),xx,n)
        do j=1,n
          do i=1,n
            do k=1,n
              y(i,j,ii)=y(i,j,ii)+conjg(xloc(k,i))*xx(k,j)
            enddo
          enddo
        enddo
      enddo
      lx=(0.d0,0.d0)
      do i_tau=0,n_tau
        do j=1,n
          de=eloc(j)-chem_pot
	    gt=green0(de,tau_mesh(i_tau))
	    do i=1,n
	      de1=eloc(j)-eloc(i)
	      if(i_nu==0.and.abs(de1)<1.d-10) then
	        gt0=green0(de,betta_t)
	        ll(i,j)=gt*(betta_t*gt0+tau_mesh(i_tau))/2.d0
	      else
	        ll(i,j)=-gt/dcmplx(de1,w_nu(i_nu))
	      endif
	    enddo
	  enddo
	  do ii=1,np
	    do j=1,n
	      do i=1,n
	        do l=1,n
	          do k=1,n
	            lx(i,j,ii,i_tau)=lx(i,j,ii,i_tau)+xloc(i,k)*y(k,l,ii)
     &	                                           *ll(k,l)
     &                                             *conjg(xloc(j,l))
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
	deallocate(ll,xx,y)
      end
      
      
      subroutine l_g1_w1_tau_x_rel(eloc,xloc,n,np,lx,i_nu,ffb)
      use atom_mod
	use heg_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: n,np,i_nu
      real*8, intent(in) :: eloc(n)
      complex*16, intent(in) :: xloc(n,n),ffb(n,n,np)
      complex*16, intent(out) :: lx(n,n,np,0:n_tau)
      integer :: i,j,ii,i_tau,l,k
      real*8 :: de,gt,green0,gt0,de1
      complex*16, allocatable :: ll(:,:),xx(:,:),y(:,:,:)
      allocate(xx(n,n))
      allocate(ll(n,n))
      allocate(y(n,n,np))
      do ii=1,np
	  call zgemm('c','n',n,n,n,(1.d0,0.d0),xloc,n,ffb(1,1,ii),n,
     &	           (0.d0,0.d0),xx,n)
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),xx,n,xloc,n,(0.d0,0.d0),
     &	           y(1,1,ii),n)
      enddo
      lx=(0.d0,0.d0)
      do i_tau=0,n_tau
        do j=1,n
          de=eloc(j)-chem_pot
	    gt=green0(de,tau_mesh(i_tau))
	    do i=1,n
	      de1=eloc(j)-eloc(i)
	      if(i_nu==0.and.abs(de1)<1.d-10) then
	        gt0=green0(de,betta_t)
	        ll(i,j)=gt*(betta_t*gt0+tau_mesh(i_tau))/2.d0
	      else
	        ll(i,j)=-gt/dcmplx(de1,w_nu(i_nu))
	      endif
	    enddo
	  enddo
	  do ii=1,np
	    do j=1,n
	      do i=1,n
	        do l=1,n
	          do k=1,n
	            lx(i,j,ii,i_tau)=lx(i,j,ii,i_tau)+xloc(i,k)*y(k,l,ii)
     &	                                           *ll(k,l)*xloc(j,l)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
	deallocate(ll,xx,y)
      end