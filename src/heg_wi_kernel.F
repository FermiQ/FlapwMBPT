      subroutine heg_wi_kernel(iq,kap,krnl)
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      implicit none
      integer, intent(in) :: iq
      real*8, intent(in) :: kap
      real*8, intent(out) :: krnl(0:n_tau)
      integer :: i_tau,ind_nu,ind_tau,i_nu
      real*8 :: vv,akap,green0_boson_tau,dd,tau,a,b,t1,t2,t3,v2
      real*8, allocatable :: tt(:,:),t(:)
      complex*16, allocatable :: ab(:)
      akap=abs(kap)
! -------------- Static part ----------------------------------
      vv=q_m_heg(iq)**2*v_q_heg(iq)
      if(abs(kap)<1.d-10) then
        do i_tau=0,n_tau
          krnl(i_tau)=-vv*(tau_mesh(i_tau)/betta_t-0.5d0)
        enddo
      else
        do i_tau=0,n_tau
          tau=tau_mesh(i_tau)
          dd=kap*(green0_boson_tau(akap,tau)-1.d0/betta_t/akap/akap)
          t1=akap*tau
          t2=akap*betta_t
	    t3=akap*(betta_t-tau)
	    a=0.d0
	    if(t1<100.d0) a=exp(-t1)
	    if(t3<100.d0) a=a-exp(-t3)
	    b=1.d0
	    if(t2<100.d0) b=b-exp(-t2)
	    dd=dd+0.5d0*a/b
          krnl(i_tau)=vv*dd
        enddo
      endif
! ---- Dynamic part --------------------------------------------
      allocate(ab(ndim3_nu))
      allocate(tt(2,ndim3_tau))
      ab=(0.d0,0.d0)
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(i_nu==0) cycle
        v2=w_wi_q_nu_heg(iq,i_nu)/(w_nu(i_nu)**2+kap*kap)
        ab(ind_nu)=dcmplx(kap*v2,w_nu(i_nu)*v2)
      enddo
c      call from_nu_to_tau_c(ab,tt,1)
      allocate(t(0:n_tau))
      t=0.d0
      do ind_tau=1,ndim3_tau
        i_tau=me_t*ndim3_tau+ind_tau-1
        t(i_tau)=tt(1,ind_tau)
        t(n_tau-i_tau)=tt(2,ind_tau)
      enddo
	if(nproc_t/=1) call dgop(t,n_tau+1,'  +',comm_t)
      krnl=krnl+q_m_heg(iq)**2*t
      deallocate(ab,tt,t)
      end