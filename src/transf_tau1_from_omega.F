      subroutine transf_tau1_from_omega
c     Creates transformation matrices for OMEGA;TAU ---> TAU1;TAU
c       They correspond to the transforms :
c     F(t-t';t) = 2/b Sum_w>0 
c           { cos(wt') Re(A(w;t) - sin(wt') Im(A(w;t) 
c           + cos(w(t-t')) Re(B(w;t) + sin(w(t-t')) Im(B(w;t) }
	use manager_mod
	use parallel_mod
	use units_mod
	use vertex_mod
      implicit none
	integer :: ind_tau,i_omega,i_tau,n,ic_omega,i1_tau,it
	real*8 :: c,c1,cs,sn,w0,wm,a,b,dd,w,x
	real*8, allocatable :: cheb(:),tmp1(:,:,:)
c ----------------------------------------------------------------
      tau1_from_omega=0.d0
      w0=sqrt(w_omega(n_omega_exa))
      wm=sqrt(omega_max)
      a=2.d0*w0*wm
      b=w0+wm
      c=wm-w0
	allocate(cheb(0:n_omega_asy-2))
	allocate(tmp1(0:n_omega_asy-2,0:n_tau,4))
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
        do it=1,2
          if(it==2) i_tau=n_tau-i_tau
		  tmp1=0.d0
	    dd=2.d0*pi/betta_t
	    w=w_omega(n_omega_exa)
	    do i_omega=n_omega_exa+1,n_omega_maz
	      w=w+dd
            x=(b-a/sqrt(w))/c
	      call dcheb_recurr(x,cheb,n_omega_asy-2)
            do i1_tau=0,n_tau1(i_tau)
	        c1=w*tau1_mesh(i1_tau,i_tau)
	        cs=cos(c1)
	        sn=sin(c1)
	        do n=0,n_omega_asy-2
	          tmp1(n,i1_tau,1)=tmp1(n,i1_tau,1)+cs*cheb(n)
	          tmp1(n,i1_tau,2)=tmp1(n,i1_tau,2)+sn*cheb(n)
	        enddo
	        c1=w*(tau_mesh(i_tau)-tau1_mesh(i1_tau,i_tau))
	        cs=cos(c1)
	        sn=sin(c1)
	        do n=0,n_omega_asy-2
	          tmp1(n,i1_tau,3)=tmp1(n,i1_tau,3)+cs*cheb(n)
	          tmp1(n,i1_tau,4)=tmp1(n,i1_tau,4)+sn*cheb(n)
	        enddo
	      enddo
	    enddo
          do i_omega=0,n_omega
	      w=w_omega(i_omega)
            if(i_omega<=n_omega_exa) then
              do i1_tau=0,n_tau1(i_tau)
	          c1=w*tau1_mesh(i1_tau,i_tau)
	          tau1_from_omega(i_omega,i1_tau,1,1,it,ind_tau)=cos(c1)
	          tau1_from_omega(i_omega,i1_tau,2,1,it,ind_tau)=sin(c1)
	          c1=w*(tau_mesh(i_tau)-tau1_mesh(i1_tau,i_tau))
	          tau1_from_omega(i_omega,i1_tau,1,2,it,ind_tau)=cos(c1)
	          tau1_from_omega(i_omega,i1_tau,2,2,it,ind_tau)=sin(c1)
	        enddo
            else if(i_omega<=n_omega-1) then
              ic_omega=i_omega-n_omega_exa
              do i1_tau=0,n_tau1(i_tau) 
	          tau1_from_omega(i_omega,i1_tau,1,1,it,ind_tau)=
     &	        dot_product(tmp1(:,i1_tau,1),cheb_omega(:,ic_omega))
	          tau1_from_omega(i_omega,i1_tau,2,1,it,ind_tau)=
     &	        dot_product(tmp1(:,i1_tau,2),cheb_omega(:,ic_omega))
	          tau1_from_omega(i_omega,i1_tau,1,2,it,ind_tau)=
     &	        dot_product(tmp1(:,i1_tau,3),cheb_omega(:,ic_omega))
	          tau1_from_omega(i_omega,i1_tau,2,2,it,ind_tau)=
     &	        dot_product(tmp1(:,i1_tau,4),cheb_omega(:,ic_omega))
	        enddo
	      endif
	    enddo
	  enddo
	enddo
      deallocate(tmp1,cheb)
c ---------- Normalization -----------------------------------------
      tau1_from_omega=tau1_from_omega*2.d0/betta_t
      end