	subroutine output_sig_vrt_bands(err_lam_tau0_nu,s_tau,lam_tau0_nu)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: err_lam_tau0_nu(0:n_nu_exa,nqdiv,2),
     &                      s_tau(nbndf,0:n_tau,npnt,nspin),
     &                      lam_tau0_nu(0:n_nu,npnt,2)
	integer :: ispin,n,i,ii,n0,m,i_len,l1,l2,i_omega,k,iq,i_nu
	complex*16, allocatable :: s_om(:,:,:,:)
	allocate(s_om(nbndf,0:n_omega,npnt,nspin))
	do ispin=1,nspin
	  do k=1,npnt
c	    call from_tau_to_omega(s_om(1,0,k,ispin),s_tau(1,0,k,ispin),
c     &	                         nbndf)
        enddo
	enddo
	n0=max(n_bnd(1,1),n_bnd(1,nspin))
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	if(maswrk) then
	  i_len=len_trim(allfile)
	  open(3,file=allfile(1:i_len)//'_Im_FULL_Vrt_Sigma_Gamma')
	  open(4,file=allfile(1:i_len)//'_Re_FULL_Vrt_Sigma_Gamma')
	endif
	do ispin=1,nspin
	  if(maswrk) write(3,*)' Ispin = ',ispin
	  if(maswrk) write(4,*)' Ispin = ',ispin
	  do i=1,n
          l1=m*(i-1)+1
          l2=min(m*i,n0)
	    if(maswrk) then
	      write(3,*)' Bands from ',l1,' to ',l2
	      write(4,*)' Bands from ',l1,' to ',l2
	      do i_omega=0,n_omega
	        write(3,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*evolt/2.d0,
     &		  (imag(s_om(ii,i_omega,1,ispin))*evolt/2.d0,ii=l1,l2)
	        write(4,'(f11.4,5(1x,f12.7))')w_omega(i_omega)*evolt/2.d0,
     &		     (real(s_om(ii,i_omega,1,ispin))*evolt/2.d0,ii=l1,l2)
	      enddo
	    endif
	  enddo
	enddo   !! over ispin
	if(maswrk) close(3)
	if(maswrk) close(4)
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_FULL_Sigma_Gamma_Tau_0')
	  do i=1,nbndf
	    write(3,'(i4,2(1x,f12.7))')i,(s_tau(i,0,1,ispin)*evolt/2.d0,
     &	                                ispin=1,nspin)
	  enddo
	  close(3)
	endif
	deallocate(s_om)
c -------------------------------------------------------------------	
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_ERR_LAM_tau0_nu')
	  do i_nu=0,n_nu_exa
	    write(3,'(e12.5,3(1x,e9.2,a1,e9.2))')w_nu(i_nu)*evolt/2.d0,
     &	    (err_lam_tau0_nu(i_nu,iq,1),'/',
     &       err_lam_tau0_nu(i_nu,iq,2),iq=1,min(3,nqdiv))
        enddo
	  close(3)
	endif
c -------------------------------------------------------------------	
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_LAM_tau0_nu_stat')
	  open(4,file=allfile(1:i_len)//'_LAM_tau0_nu_dyn')
	  do i_nu=0,n_nu
	    write(3,'(e12.5,6(1x,e9.2))')w_nu(i_nu)*evolt/2.d0,
     &	    (lam_tau0_nu(i_nu,iq,1),iq=1,min(6,npnt))
	    write(4,'(e12.5,6(1x,e9.2))')w_nu(i_nu)*evolt/2.d0,
     &	    (lam_tau0_nu(i_nu,iq,2),iq=1,min(6,npnt))
        enddo
	  close(3)
	  close(4)
	endif
      end