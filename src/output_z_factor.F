      subroutine output_z_factor(zd)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      real*8, intent(in) :: zd(nbndf,ndim3_k(me_k+1),nspin)
	  integer :: k,i_len,m,n,ispin,l1,l2,m0,ind_k,i,nb,ib,n0
	  real*8, allocatable :: z_out(:,:)
	  n0=nbndf
	  m=5
	  n=npnt/m
	  if(n*m<npnt) n=n+1
	  allocate(z_out(n0,npnt))
	  if(maswrk) then
	    i_len=len_trim(allfile)
	    if(ubi==' gw') then
	      open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.gw')
	    else if(ubi==' qp') then
	      open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.qp')
	    else if(ubi=='psi') then
	      open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.psi')
	    else if(ubi=='bsp') then
	      open(3,file=allfile(1:i_len)//'_Z_factor_band_k_x.bsp')
	    endif
	  endif
	  do ispin=1,nspin
	    if(maswrk) write(3,*)' Ispin = ',ispin
	    z_out=0.d0
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      nb=n_bnd(k,ispin)
	      do i=1,nb
	        z_out(i,k)=zd(i,ind_k,ispin)
	      enddo
	    enddo
c --------- Collect from all processes  ------------------------------------------
	    if(nproc_k/=1) call dgop(z_out,n0*npnt,'  +',comm_k)
	    do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,npnt)
	      if(maswrk) write(3,*)' K-points from ',l1,' to ',l2
	      m0=l2-l1+1
	      if(maswrk) then
	        do ib=1,n0
	          write(3,'(i5,5(1x,f12.7))')ib,(z_out(ib,k),k=1,m0)
	        enddo
	      endif
	    enddo
	  enddo
	  if(maswrk) then
	    close(3)
	  endif
	  deallocate(z_out)
      end
