      subroutine from_nu_omega_to_tau_omega_2(i_omega,fnu,ftau,n,fnux)
c     Makes the 1/nu**2 correction for particle-particle case      
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: i_omega,n
      real*8, intent(in) :: fnux(n)
      complex*16, intent(inout) :: fnu(n,n_nu1_max)
	complex*16, intent(out) :: ftau(n,2,ndim3_tau,2)
      integer :: nn,n0,i_tau,it,ind_tau,i_nu
      real*8 :: phase,ff,vv,ww,ar,ai
      complex*16 :: cc
      n0=nu_adapt_index_zero(i_omega)  ! index for nu = 0
      ww=w_omega(i_omega)
c ----- Subtract nu/=0 terms ----------------------------------
      do i_nu=1,num_nu1_adapt(i_omega)
        if(i_nu==n0) cycle
        vv=w_nu_adapt_om(i_nu,i_omega)
        ff=1.d0/(ww-vv)/vv
        fnu(:,i_nu)=fnu(:,i_nu)-fnux*ff
      enddo
c -------------------- A - component ( v <= 0 ) ------------------
      nn=n0	
	call zgemm('n','n',n,2*ndim3_tau,nn,(1.d0,0.d0),fnu,n,
     &	         tau_from_nu_omega(1,1,1,i_omega),n_nu1_max,(0.d0,0.d0),
     &           ftau(1,1,1,1),n)
c -------------------- B - component ( v >= 0 ) ------------------
      nn=num_nu1_adapt(i_omega)-n0+1
	call zgemm('n','n',n,2*ndim3_tau,nn,(1.d0,0.d0),fnu(1,n0),n,
     &	         tau_from_nu_omega(n0,1,1,i_omega),n_nu1_max,
     &           (0.d0,0.d0),ftau(1,1,1,2),n)
	do ind_tau=1,ndim3_tau
	  do it=1,2
	    i_tau=me_t*ndim3_tau+ind_tau-1
	    if(it==2) i_tau=n_tau-i_tau
          phase=-w_omega(i_omega)*tau_mesh(i_tau)
          cc=dcmplx(cos(phase),sin(phase))
          ftau(:,it,ind_tau,2)=cc*ftau(:,it,ind_tau,2)
        enddo
      enddo
c ------- Return the leading terms ----------------------
	do ind_tau=1,ndim3_tau
	  do it=1,2
	    i_tau=me_t*ndim3_tau+ind_tau-1
	    if(it==2) i_tau=n_tau-i_tau
	    if(i_tau==0) cycle
	    if(i_tau==n_tau) cycle
	    ar=-1.d0/betta_t/ww/ww
	    ai=(0.5d0-tau_mesh(i_tau)/betta_t)/ww
          ftau(:,it,ind_tau,1)=ftau(:,it,ind_tau,1)+fnux*dcmplx(ar,ai)
	    ar=0.d0
	    ai=-1.d0/2.d0/ww
          ftau(:,it,ind_tau,2)=ftau(:,it,ind_tau,2)+fnux*dcmplx(ar,ai)
        enddo
      enddo      
      end