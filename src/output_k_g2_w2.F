      subroutine output_k_g2_w2(key,k_pw,i_nu,n,nq)
c     key = 1  for K0
c     key = 2  for K      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: key,i_nu,n,nq
	complex*16, intent(in) :: k_pw(n,n,0:n_tau,2,nq,nspin)
	integer :: i_len,i_tau,np,k,ispin
	real*8 :: phase
	complex*16 :: cc,p(5)
	i_len=len_trim(allfile)
	np=min(5,nq)
	if(maswrk) then
	  if(key==1) then
	    open(3,file=allfile(1:i_len)//'_K_0_vertex_Re')
	    open(4,file=allfile(1:i_len)//'_K_0_vertex_Im')
	    open(5,file=allfile(1:i_len)//'_K_0_AK_BK_comp')
	  else if(key==2) then
	    open(3,file=allfile(1:i_len)//'_K_vertex_Re')
	    open(4,file=allfile(1:i_len)//'_K_vertex_Im')
	    open(5,file=allfile(1:i_len)//'_K_AK_BK_comp')
	  endif
	  do ispin=1,nspin
	    write(3,*)' Ispin = ',ispin
	    write(4,*)' Ispin = ',ispin
	    write(5,*)' Ispin = ',ispin
	    do i_tau=0,n_tau
            phase=w_nu(i_nu)*tau_mesh(i_tau)
            cc=dcmplx(cos(phase),sin(phase))
            do k=1,np
              p(k)=k_pw(1,1,i_tau,1,k,ispin)
     &            +cc*k_pw(1,1,i_tau,2,k,ispin)
            enddo
	      write(3,'(f11.5,5(1x,e12.5))')tau_mesh(i_tau),
     &	                                  (real(p(k)),k=1,np)
	      write(4,'(f11.5,5(1x,e12.5))')tau_mesh(i_tau),
     &	                                  (imag(p(k)),k=1,np)
	      write(5,'(f11.5,6(1x,e10.3))')tau_mesh(i_tau),
     &	                                  real(p(1)),
     &                         real(k_pw(1,1,i_tau,1,1,ispin)),
     &                         real(k_pw(1,1,i_tau,2,1,ispin)),
     &	                                  imag(p(1)),
     &                         imag(k_pw(1,1,i_tau,1,1,ispin)),
     &                         imag(k_pw(1,1,i_tau,2,1,ispin))
          enddo
	  enddo
	  close(3)
	  close(4)
	  close(5)
	endif
      end