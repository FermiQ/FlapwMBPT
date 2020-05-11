	subroutine heg_w_from_resp
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      use vertex_mod
      use units_mod
	implicit none
	integer :: i,i_nu,ind_tau,i_tau,ind_nu,iq
	real*8 :: e,vp,p,pi2a,q,x,q2
	real*8, allocatable :: e1(:,:),wqtau(:,:),wqtauc(:,:)
      pi2a=(pi+pi)/par
	allocate(e1(0:n_q_heg,0:n_nu))
	e1=0.d0
	do i=n_q_heg,0,-1
	  if(i==0) then
	    call heg1_e1_q0_fit(e1)
	    p=1.d0/resp_q_nu_heg(i,0)+v_q_heg(0)
	    p=1.d0/p
	    w_q_nu_heg(0,0)=-v_q_heg(0)
	    if(iter_psi+iter_bsp/=0) w_q_nu_heg_c(0,1)=-v_q_heg_c(1)
	    p_q_nu_heg(0,0)=p
          edd_q_nu_heg(0,0)=1.d3
	    do i_nu=1,n_nu
	      e=1.d0/e1(0,i_nu)
            edd_q_nu_heg(0,i_nu)=e
	      w_q_nu_heg(0,i_nu)=v_q_heg(0)*e1(0,i_nu)-v_q_heg(0)
	      p_q_nu_heg(0,i_nu)=resp_q_nu_heg(0,i_nu)/e1(0,i_nu)
	      if(iter_psi+iter_bsp/=0) w_q_nu_heg_c(i_nu,1)=
     &                     v_q_heg_c(1)*e1(0,i_nu)-v_q_heg_c(1)
	    enddo
	  else
	    do i_nu=0,n_nu
	      vp=v_q_heg(i)*resp_q_nu_heg(i,i_nu)
	      e1(i,i_nu)=1.d0+vp
            edd_q_nu_heg(i,i_nu)=1.d0/e1(i,i_nu)
	      w_q_nu_heg(i,i_nu)=vp*v_q_heg(i)
	      p_q_nu_heg(i,i_nu)=resp_q_nu_heg(i,i_nu)/e1(i,i_nu)
	    enddo
	  endif
	enddo
	do i=0,n_q_heg
	  do i_nu=0,n_nu
	    e1_heg(i,i_nu)=e1(i,i_nu)
	  enddo
	enddo
c ----------- W_Q_Tau ----------------------------------------  
      do i=0,n_q_heg
        call nu_to_tau_spl(w_q_nu_heg(i,:),w_q_tau_heg(i,:),1)
        call nu_to_tau_spl(resp_q_nu_heg(i,:),resp_q_tau_heg(i,:),1)
      enddo
      w_q_tau_heg=w_q_tau_heg/2
      resp_q_tau_heg=resp_q_tau_heg/2
c -------- Splines for W(q,nu) -------------------------------- 
      do i_nu=0,n_nu
        call spline_inhmg(q_m_heg,w_q_nu_heg(0,i_nu),
     &                    w_q_nu_spl(1,1,i_nu),w_q_nu_spl(1,2,i_nu),
     &                    w_q_nu_spl(1,3,i_nu),w_q_nu_spl(1,4,i_nu),
     &                    n_q_heg,0,0.d0,0.d0)
      enddo
      if(iter_psi+iter_bsp/=0) then
        allocate(wqtau(0:n_tau/2,nstar_c))
        allocate(wqtauc(0:n_tau/2,nstar_c))    
        do i=1,nstar_c
          if(i/=1) then
            q=pi2a*qabs_star_c(i)
            q2=q*q
            x=sqrt(q)
            iq=x/step_q_x
            if(iq*step_q_x<x) iq=iq+1
            if(iq==0) iq=1
            q=q-q_m_heg(iq-1)
            do i_nu=0,n_nu
              w_q_nu_heg_c(i_nu,i)=w_q_nu_spl(iq,1,i_nu)+q*(
     &                             w_q_nu_spl(iq,2,i_nu)+q*(
     &                             w_q_nu_spl(iq,3,i_nu)+q*
     &                             w_q_nu_spl(iq,4,i_nu)))
            enddo
          endif
          call nu_to_tau_spl(w_q_nu_heg_c(0,i),wqtau(0,i),1)
          call nu_to_tau_spl(w_q_nu_heg_c(0,i),wqtauc(0,i),1)
        enddo
        wqtau=wqtau/2
        wqtauc=wqtauc/2
        w_q_tau_heg_c=wqtau
        w_r_tau_heg_c=0.d0
	  do ind_tau=1,ndim3_tau
	    i_tau=ndim3_tau*me_t+ind_tau-1
	    call from_q_to_rr_heg_red(nstar_c,wqtauc(i_tau,:),
     &	                          w_r_tau_heg_c(1,i_tau))
	  enddo
	  if(nproc_t/=1) call dgop(w_r_tau_heg_c,
     &                             nrr_div_c*(n_tau/2+1),'  +',
     &                             comm_t)
        deallocate(wqtau,wqtauc)
	endif
c ----------- W_R_Tau ----------------------------------------
      w_r_tau_heg=0.d0
      resp_r_tau_heg=0.d0            
	do ind_tau=1,ndim3_tau
	  i_tau=me_t*ndim3_tau+ind_tau-1
	  call from_q_to_r_heg1(w_q_tau_heg(0,i_tau),
     &	                    w_r_tau_heg(0,i_tau))
        resp_q_tau_heg(:,i_tau)=resp_q_tau_heg(:,i_tau)+ro_el_gas
	  call from_q_to_r_heg1(resp_q_tau_heg(0,i_tau),
     &	                    resp_r_tau_heg(0,i_tau))
        resp_q_tau_heg(:,i_tau)=resp_q_tau_heg(:,i_tau)-ro_el_gas
	enddo
	if(nproc_t/=1) then
	  call dgop(w_r_tau_heg,(n_r_heg+1)*(n_tau/2+1),'  +',comm_t)
	  call dgop(resp_r_tau_heg,(n_r_heg+1)*(n_tau/2+1),'  +',comm_t)
	endif
c ----------- W_R_Nu ----------------------------------------
      w_r_nu_heg=0.d0      
	do ind_nu=1,ndim3_nu
	  i_nu=me_t*ndim3_nu+ind_nu-1
	  call from_q_to_r_heg1(w_q_nu_heg(0,i_nu),
     &	                    w_r_nu_heg(0,i_nu))
	enddo
	if(nproc_t/=1) call dgop(w_r_nu_heg,(n_r_heg+1)*(n_nu+1),'  +',
     &                           comm_t)
c -----------------------------------------------------------------     
	call output_heg1_v_r(1)
	call output_heg1_g_pair
	deallocate(e1)
      end
      
      
      subroutine heg1_w   
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      use vertex_mod
	implicit none
	integer :: i,i_nu,ind_tau,i_tau,ind_nu,iq
	real*8 :: e,vp,q,pi2a,x,q2,e11,roeff,a,kp2,v,q4
	real*8, allocatable :: e1(:,:),wqtau(:,:),wqtauc(:,:)
	pi2a=(pi+pi)/par
	allocate(e1(0:n_q_heg,0:n_nu))
	e1=0.d0
c ---------------------------------------------------------------      
      do i_nu=n_nu,0,-1
        v=w_nu(i_nu)
	  do i=n_q_heg,0,-1
	    q=q_m_heg(i)
          q4=q**4
          if(q>=q_small_heg) then
	      vp=v_q_heg(i)*p_q_nu_heg(i,i_nu)
	      e=1.d0-vp
            edd_q_nu_heg(i,i_nu)=e
            e1(i,i_nu)=1.d0/e
            e11=e1(i,i_nu)-1.d0
            if(q_m_heg(i-1)<q_small_heg) then
              call ro_eff_finder(e11,q,v,roeff)
            endif
	      resp_q_nu_heg(i,i_nu)=p_q_nu_heg(i,i_nu)*e1(i,i_nu)
          else
            a=16.d0*pi*roeff
            kp2=q4+a
            e11=-a/(v*v+kp2)
            e1(i,i_nu)=e11+1.d0
            edd_q_nu_heg(i,i_nu)=1.d0/e1(i,i_nu)
            if(i_nu/=0) then
              resp_q_nu_heg(i,i_nu)=-2.d0*roeff*q*q/(v*v+kp2)
              p_q_nu_heg(i,i_nu)=-2.d0*roeff*q*q/(v*v+q4)
            endif
          endif
	    w_q_nu_heg(i,i_nu)=e11*v_q_heg(i)
        enddo
      enddo
c ---------------------------------------------------------------        
	do i=0,n_q_heg
	  do i_nu=0,n_nu
	    e1_heg(i,i_nu)=e1(i,i_nu)
	  enddo
	enddo
c ----------- W_Q_Tau ----------------------------------------
      call nu_to_tau_spl(w_q_nu_heg,w_q_tau_heg,n_q_heg+1)
      call nu_to_tau_spl(resp_q_nu_heg,resp_q_tau_heg,n_q_heg+1)
      w_q_tau_heg=w_q_tau_heg/2.d0
      resp_q_tau_heg=resp_q_tau_heg/2.d0
      do i=0,n_q_heg
        dw_q_tau_heg(i)=-w_q_nu_heg(i,n_nu)*w_nu(n_nu)**2/2
      enddo
c -------- Splines for W(q,nu) -------------------------------- 
      do i_nu=0,n_nu
        call spline_inhmg(q_m_heg,w_q_nu_heg(0,i_nu),
     &                    w_q_nu_spl(1,1,i_nu),w_q_nu_spl(1,2,i_nu),
     &                    w_q_nu_spl(1,3,i_nu),w_q_nu_spl(1,4,i_nu),
     &                    n_q_heg,0,0.d0,0.d0)
      enddo
      if(iter_psi+iter_bsp/=0) then
	  w_q_nu_heg_c(0,1)=-v_q_heg_c(1)
	  do i_nu=1,n_nu
	    w_q_nu_heg_c(i_nu,1)=v_q_heg_c(1)*e1(0,i_nu)-v_q_heg_c(1)
        enddo
        allocate(wqtau(0:n_tau/2,nstar_c))
        allocate(wqtauc(0:n_tau/2,nstar_c))    
        do i=1,nstar_c
          if(i/=1) then
            q=pi2a*qabs_star_c(i)
            q2=q*q
            x=sqrt(q)
            iq=x/step_q_x
            if(iq*step_q_x<x) iq=iq+1
            if(iq==0) iq=1
            q=q-q_m_heg(iq-1)
            do i_nu=0,n_nu
              w_q_nu_heg_c(i_nu,i)=w_q_nu_spl(iq,1,i_nu)+q*(
     &                             w_q_nu_spl(iq,2,i_nu)+q*(
     &                             w_q_nu_spl(iq,3,i_nu)+q*
     &                             w_q_nu_spl(iq,4,i_nu)))
            enddo
          endif
          call nu_to_tau_spl(w_q_nu_heg_c(0,i),wqtau(0,i),1)
          wqtau(:,i)=wqtau(:,i)/2.d0
          call nu_to_tau_spl(w_q_nu_heg_c(0,i),wqtauc(0,i),1)
          wqtauc(:,i)=wqtauc(:,i)/2.d0
        enddo
        w_q_tau_heg_c=wqtau
        w_r_tau_heg_c=0.d0
	  do ind_tau=1,ndim3_tau
	    i_tau=ndim3_tau*me_t+ind_tau-1
	    call from_q_to_rr_heg_red(nstar_c,wqtauc(i_tau,:),
     &	                          w_r_tau_heg_c(1,i_tau))
	  enddo
	  if(nproc_t/=1) call dgop(w_r_tau_heg_c,
     &                             nrr_div_c*(n_tau/2+1),'  +',
     &                             comm_t)
        deallocate(wqtau,wqtauc)
	endif
c ----------- W_R_Tau ----------------------------------------
      w_r_tau_heg=0.d0
      resp_r_tau_heg=0.d0            
	do ind_tau=1,ndim3_tau
	  i_tau=me_t*ndim3_tau+ind_tau-1
	  call from_q_to_r_heg1(w_q_tau_heg(0,i_tau),
     &	                    w_r_tau_heg(0,i_tau))
        resp_q_tau_heg(:,i_tau)=resp_q_tau_heg(:,i_tau)+ro_el_gas
	  call from_q_to_r_heg1(resp_q_tau_heg(0,i_tau),
     &	                    resp_r_tau_heg(0,i_tau))
        resp_q_tau_heg(:,i_tau)=resp_q_tau_heg(:,i_tau)-ro_el_gas
	enddo
	if(nproc_t/=1) then
	  call dgop(w_r_tau_heg,(n_r_heg+1)*(n_tau/2+1),'  +',comm_t)
	  call dgop(resp_r_tau_heg,(n_r_heg+1)*(n_tau/2+1),'  +',comm_t)
	endif
c ----------- W_R_Nu ----------------------------------------
      w_r_nu_heg=0.d0      
	do ind_nu=1,ndim3_nu
	  i_nu=me_t*ndim3_nu+ind_nu-1
	  call from_q_to_r_heg1(w_q_nu_heg(0,i_nu),
     &	                    w_r_nu_heg(0,i_nu))
	enddo
	if(nproc_t/=1) call dgop(w_r_nu_heg,(n_r_heg+1)*(n_nu+1),'  +',
     &                           comm_t)
c -----------------------------------------------------------------     
	call output_heg1_v_r(1)
	call output_heg1_g_pair
	deallocate(e1)
	end
      
      
      subroutine ro_eff_finder(e11,q,v,roeff)
c ------- MESH for Nu ------------------------------------------------
	use manager_mod
	use models_mod
	use parallel_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: e11,q,v
      real*8, intent(out) :: roeff
	real*8 :: q1,q2,q3,f1,f2,f3,kp2,a,q4
      q4=q**4
      q1=0.d0
      q2=1.d5
      a=16.d0*pi*q1
      kp2=q4+a
      f1=e11+a/(v*v+kp2)
      a=16.d0*pi*q2
      kp2=q4+a
      f2=e11+a/(v*v+kp2)
1     q3=(q1+q2)/2
      a=16.d0*pi*q3
      kp2=q4+a
      f3=e11+a/(v*v+kp2)
      if(f3<0.d0) then
        q1=q3
        f1=f3
      else
        q2=q3
        f2=f3
      endif
      if(f2-f1>1.d-12) goto 1
      roeff=(q1+q2)/2
      end