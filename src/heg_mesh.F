      subroutine heg_mesh
	use heg_mod
	use manager_mod
	use parallel_mod
      use solid_mod
	use units_mod
      implicit none
	integer :: i,ik,ir,iq
	real*8 :: xi,xmax,k,r,q,kmi,kma,pi2a
      pi2a=(pi+pi)/par
c --------- K-mesh -----------------------------------------------------
      allocate(k_m_heg(0:n_k_heg))
      allocate(k_w_heg(0:n_k_heg))
      x_fermi=asinh(k_fermi/alf_k_heg)
      x_k_infinum=x_fermi+asinh((k_infinum-k_fermi)/alf_k_heg)
      step_k_x=x_k_infinum/n_k_heg
      n_k_fermi=x_fermi/step_k_x
      step_k_x=x_fermi/n_k_fermi
      k_m_heg(0)=0.d0
      k_w_heg(0)=alf_k_heg*cosh(x_fermi)
      do i=1,n_k_heg
        k_m_heg(i)=k_fermi+alf_k_heg*sinh(step_k_x*i-x_fermi)
        k_w_heg(i)=alf_k_heg*cosh(step_k_x*i-x_fermi)
        if(mod(i,2)==0.and.k_m_heg(i)<=k_infinum) ik=i
      enddo
      n_k_heg=ik
      x_k_infinum=step_k_x*n_k_heg
      k_infinum=k_fermi+alf_k_heg*sinh(x_k_infinum-x_fermi)
      if(iter_psi+iter_bsp/=0) then
        if(k_infinum<=pi2a*maxval(qabs_star_c(1:nstar_c))) then
          if(maswrk) write(iun,*)' K_INFINUM < Q_MAX_CUBE '
          call ending
        endif
      endif
      kmi=1.d5
      kma=0.d0
      do i=1,n_k_heg
          kmi=min(kmi,k_m_heg(i)-k_m_heg(i-1))
          kma=max(kma,k_m_heg(i)-k_m_heg(i-1))
      enddo
      if(maswrk) write(iun,*)' Smallest K-mesh interval: ',kmi
      if(maswrk) write(iun,*)' Largest K-mesh interval: ',kma
c ---------------- 2 points ----------------------------------------
c      k_w_heg(0)=k_w_heg(0)*step_k_x/2.d0
c      k_w_heg(n_k_heg)=k_w_heg(n_k_heg)*step_k_x/2.d0
c ---------------- 3 points ----------------------------------------
      k_w_heg(0)=k_w_heg(0)*step_k_x/3.d0
      k_w_heg(n_k_heg)=k_w_heg(n_k_heg)*step_k_x/3.d0
c ---------------- 7 points ----------------------------------------
c      k_w_heg(0)=k_w_heg(0)*41.d0/140.d0*step_k_x
c      k_w_heg(n_k_heg)=k_w_heg(n_k_heg)*41.d0/140.d0*step_k_x
      do i=1,n_k_heg-1
c ---------------- 2 points ----------------------------------------
c        k_w_heg(i)=k_w_heg(i)*step_k_x
c ---------------- 3 points ----------------------------------------
       if(mod(i,2)==1) k_w_heg(i)=k_w_heg(i)*4.d0/3.d0*step_k_x
       if(mod(i,2)==0) k_w_heg(i)=k_w_heg(i)*2.d0/3.d0*step_k_x
c ---------------- 7 points ----------------------------------------
c        if(mod(i,6)==1) k_w_heg(i)=k_w_heg(i)*216.d0/140.d0*step_k_x
c        if(mod(i,6)==2) k_w_heg(i)=k_w_heg(i)*27.d0/140.d0*step_k_x
c        if(mod(i,6)==3) k_w_heg(i)=k_w_heg(i)*272.d0/140.d0*step_k_x
c        if(mod(i,6)==4) k_w_heg(i)=k_w_heg(i)*27.d0/140.d0*step_k_x
c        if(mod(i,6)==5) k_w_heg(i)=k_w_heg(i)*216.d0/140.d0*step_k_x
c        if(mod(i,6)==0) k_w_heg(i)=k_w_heg(i)*82.d0/140.d0*step_k_x
      enddo
c --------- P-mesh -----------------------------------------------------
      allocate(p_m_heg(0:n_p_heg))
      allocate(p_w_heg(0:n_p_heg))
      p_infinum=2.d0*q_infinum
      x_p_infinum=x_fermi+asinh((p_infinum-k_fermi)/alf_k_heg)
      step_p_x=x_p_infinum/n_p_heg
      n_p_fermi=x_fermi/step_p_x
      step_p_x=x_fermi/n_p_fermi
      x_p_infinum=step_p_x*n_p_heg
      p_infinum=k_fermi+alf_k_heg*sinh(x_p_infinum-x_fermi)
      p_m_heg(0)=0.d0
      p_w_heg(0)=alf_k_heg*cosh(x_fermi)
      do i=1,n_p_heg
        p_m_heg(i)=k_fermi+alf_k_heg*sinh(step_p_x*i-x_fermi)
        p_w_heg(i)=alf_k_heg*cosh(step_p_x*i-x_fermi)
      enddo
c ---------------- 2 points ----------------------------------------
c      p_w_heg(0)=p_w_heg(0)*step_p_x/2.d0
c      p_w_heg(n_p_heg)=p_w_heg(n_p_heg)*step_p_x/2.d0
c ---------------- 3 points ----------------------------------------
       p_w_heg(0)=p_w_heg(0)*step_p_x/3.d0
       p_w_heg(n_p_heg)=p_w_heg(n_p_heg)*step_p_x/3.d0
c ---------------- 7 points ----------------------------------------
c      k_w_heg(0)=k_w_heg(0)*41.d0/140.d0*step_k_x
c      k_w_heg(n_k_heg)=k_w_heg(n_k_heg)*41.d0/140.d0*step_k_x
      do i=1,n_p_heg-1
c ---------------- 2 points ----------------------------------------
c        p_w_heg(i)=p_w_heg(i)*step_p_x
c ---------------- 3 points ----------------------------------------
        if(mod(i,2)==1) p_w_heg(i)=p_w_heg(i)*4.d0/3.d0*step_p_x
        if(mod(i,2)==0) p_w_heg(i)=p_w_heg(i)*2.d0/3.d0*step_p_x
c ---------------- 7 points ----------------------------------------
c        if(mod(i,6)==1) k_w_heg(i)=k_w_heg(i)*216.d0/140.d0*step_k_x
c        if(mod(i,6)==2) k_w_heg(i)=k_w_heg(i)*27.d0/140.d0*step_k_x
c        if(mod(i,6)==3) k_w_heg(i)=k_w_heg(i)*272.d0/140.d0*step_k_x
c        if(mod(i,6)==4) k_w_heg(i)=k_w_heg(i)*27.d0/140.d0*step_k_x
c        if(mod(i,6)==5) k_w_heg(i)=k_w_heg(i)*216.d0/140.d0*step_k_x
c        if(mod(i,6)==0) k_w_heg(i)=k_w_heg(i)*82.d0/140.d0*step_k_x
      enddo
c ------------- Q mesh ---------------------------------------
      allocate(q_m_heg(0:n_q_heg))
      allocate(q_w_heg(0:n_q_heg))
	xmax=sqrt(q_infinum)
	step_q_x=xmax/n_q_heg
	do i=0,n_q_heg
        xi=step_q_x*i
	  q_m_heg(i)=xi*xi
	  q_w_heg(i)=2.d0*xi
	enddo
c ---------------- 2 points ----------------------------------------
      q_w_heg(0)=q_w_heg(0)*step_q_x/2.d0
      q_w_heg(n_q_heg)=q_w_heg(n_q_heg)*step_q_x/2.d0
      do i=1,n_q_heg-1
c ---------------- 2 points ----------------------------------------
        q_w_heg(i)=q_w_heg(i)*step_q_x
      enddo
c --------- R-mesh -----------------------------------------------------
      allocate(r_m_heg(0:n_r_heg))
      allocate(r_w_heg(0:n_r_heg))
	xmax=r_infinum**(1.d0/3.d0)
	step_r_x=xmax/n_r_heg
	do i=0,n_r_heg
        xi=step_r_x*i
	  r_m_heg(i)=xi**3
	  r_w_heg(i)=3.d0*xi*xi*step_r_x/3.d0
	enddo
      do i=1,n_r_heg-1
        if(mod(i,2)==0) r_w_heg(i)=r_w_heg(i)*2
        if(mod(i,2)/=0) r_w_heg(i)=r_w_heg(i)*4
      enddo
c --------- X-mesh -----------------------------------------------------
      allocate(x_m_heg(0:n_x_heg))
      allocate(x_w_heg(0:n_x_heg))
	xmax=2.d0
	step_x_x=xmax/n_x_heg
	do i=0,n_x_heg
        xi=step_x_x*i
	  x_m_heg(i)=xi-1.d0
	  x_w_heg(i)=0.5d0*step_x_x
	enddo
      do i=1,n_x_heg-1
        x_w_heg(i)=x_w_heg(i)*2
      enddo
c ----------------------------------------------------------      
      allocate(k_to_r(0:n_r_heg,0:n_k_heg))     
      allocate(r_to_k(0:n_k_heg,0:n_r_heg))      
      allocate(q_to_r(0:n_r_heg,0:n_q_heg))     
      allocate(r_to_q(0:n_q_heg,0:n_r_heg))
c ---------- K ----> R transform ----------------------------
      do ik=0,n_k_heg
        k=k_m_heg(ik)
        xi=k*k_w_heg(ik)/2.d0/pi/pi
        k_to_r(0,ik)=xi*k
        do ir=1,n_r_heg
          r=r_m_heg(ir)
          k_to_r(ir,ik)=xi*sin(k*r)/r
        enddo
      enddo
c ---------- Q ----> R transform ----------------------------
      do iq=0,n_q_heg
        q=q_m_heg(iq)
        xi=q*q_w_heg(iq)/2.d0/pi/pi
        q_to_r(0,iq)=xi*q
        do ir=1,n_r_heg
          r=r_m_heg(ir)
          q_to_r(ir,iq)=xi*sin(q*r)/r
        enddo
      enddo
c ---------- R ----> K transform ----------------------------
      do ir=0,n_r_heg
        r=r_m_heg(ir)
        xi=4.d0*pi*r*r_w_heg(ir)
        r_to_k(0,ir)=xi*r
        do ik=1,n_k_heg
          k=k_m_heg(ik)
          r_to_k(ik,ir)=xi*sin(k*r)/k
        enddo
      enddo
c ---------- R ----> Q transform ----------------------------
      do ir=0,n_r_heg
        r=r_m_heg(ir)
        xi=4.d0*pi*r*r_w_heg(ir)
        r_to_q(0,ir)=xi*r
        do iq=1,n_q_heg
          q=q_m_heg(iq)
          r_to_q(iq,ir)=xi*sin(q*r)/q
        enddo
      enddo
      end