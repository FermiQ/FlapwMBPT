      subroutine g_x_atom(ispin,x)
c ------ Exchange part of Green's Function is calculated here --------
c	x - Chemical Potential
      use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin
	real*8, intent(in) :: x
	integer :: ind_omega,i_omega,ind_tau,it,i_tau,li,l,ie
	real*8 :: de,f,green0
	complex*16 :: cf
	gx_omega_atom(:,:,:,ispin)=(0.d0,0.d0)
	gx_tau_atom(:,:,:,:,ispin)=0.d0
      do li=1,nrel*lmb(1)+1
	  if(irel.le.1) then
	    l=li-1
	  else
	    l=li/2
	  endif
        do ie=1,ntle(l,1)
          de=x-ex_atom(ie,li,ispin)
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
	      cf=dcmplx(de,w_omega(i_omega))
	      cf=(1.d0,0.d0)/cf
	      gx_omega_atom(ie,li,ind_omega,ispin)=cf
          enddo
          do ind_tau=1,ndim3_tau
            do it=1,2
              i_tau=me_t*ndim3_tau+ind_tau-1
              if(it==2) i_tau=n_tau-i_tau
	        f=green0(-de,tau_mesh(i_tau))
	        gx_tau_atom(ie,li,it,ind_tau,ispin)=f
            enddo
          enddo
        enddo
	enddo
      end