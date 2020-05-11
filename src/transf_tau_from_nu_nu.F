      subroutine transf_tau_from_nu_nu(key,ii,ref,n_ex,n_as,nu,tran)
c     key = -1 asymptotic is from the left
c     key =  1 asymptotic is from the right
c     ------------------------------------------------------       
c     Creates transformation matrices for NU' ---> TAU
c       NU'-mesh is adaptive for every NU ---------------
c       They correspond to the transform :
c                F(t;v) = 1/b Sum_v'  F(v')[v] * exp(iv't)
      use manager_mod
      use units_mod
      use vertex_mod
      implicit none
	integer, intent(in) :: key,ii,n_ex,n_as,nu(0:n_nu2_max)
	real*8, intent(in) :: ref
      complex*16, intent(out) :: tran(0:n_nu2_max,0:n_tau)
	integer :: i_tau,i_nu,i,i1,i2,i3
	real*8 :: bt,pi2_bet,vt,c1,c2,c3,v,nu1,nu2,nu3,dt,a1,a2,a3,b1,b2,
     &          b3,x2,x4,x6,x3,x5,x7,cs,v2,v4,v6,sn,v3,v5,v7
      complex*16 :: bti
      integer, allocatable :: mesh(:)
      allocate(mesh(0:n_ex+n_as))
      do i=0,n_ex+n_as
        if(key==-1) i1=ii+n_ex+n_as-i
        if(key==1) i1=ii+i
        mesh(i)=nu(i1)
      enddo
c -----------------------------------------------------------      
	bt=1.d0/betta_t
	bti=bt*(0.d0,1.d0)
	pi2_bet=2.d0*pi/betta_t
	do i_tau=0,n_tau
c ---------------------- Exact Part ------------------------------	  
	  do i_nu=0,n_ex
          if(key==-1) i1=ii+n_ex+n_as-i_nu
          if(key==1) i1=ii+i_nu
	    v=pi2_bet*mesh(i_nu)
	    vt=v*tau_mesh(i_tau)
	    tran(i1,i_tau)=bt*dcmplx(cos(vt),sin(vt))
	  enddo
c ----------------- Asymptotic Part ------------------------------	
	  do i=n_ex+1,n_ex+n_as-1,2
          if(key==-1) i2=ii+n_ex+n_as-i
          if(key==1) i2=ii+i
	    i1=i2-key
	    i3=i2+key
	    nu1=pi2_bet*mesh(i-1)
	    nu2=pi2_bet*mesh(i)
	    nu3=pi2_bet*mesh(i+1)
	    nu1=abs(ref-nu1)
	    nu2=abs(ref-nu2)
	    nu3=abs(ref-nu3)
c ---------------------- COS transformation -----------------------	    
	    dt=((1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2
     &	     +(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2
     &	     +(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2)
     &       /(nu1*nu2*nu3)**2
	    a1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**4/dt
	    a2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**4/dt
	    a3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**4/dt
	    b1=(1.d0/nu2**4-1.d0/nu3**4)/(nu2*nu3)**2/dt
	    b2=(1.d0/nu3**4-1.d0/nu1**4)/(nu1*nu3)**2/dt
	    b3=(1.d0/nu1**4-1.d0/nu2**4)/(nu1*nu2)**2/dt
	    c1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2/dt
	    c2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2/dt
	    c3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2/dt
	    x2=0.d0
	    x4=0.d0
	    x6=0.d0
	    do i_nu=mesh(i-1)+key,mesh(i+1),key
	      v=pi2_bet*i_nu
	      cs=cos(v*tau_mesh(i_tau))
	      v=abs(v-ref)
	      v2=v*v
	      v4=v2*v2
	      v6=v4*v2
	      x2=x2+cs/v2
	      x4=x4+cs/v4
	      x6=x6+cs/v6
	    enddo
	    tran(i1,i_tau)=tran(i1,i_tau)+(x2*a1+x4*b1+x6*c1)*bt
	    tran(i2,i_tau)=tran(i2,i_tau)+(x2*a2+x4*b2+x6*c2)*bt
	    tran(i3,i_tau)=tran(i3,i_tau)+(x2*a3+x4*b3+x6*c3)*bt
c ---------------------- SIN transformation -----------------------	    
	    dt=((1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**2
     &	     +(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**2
     &	     +(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**2)
     &       /(nu1*nu2*nu3)**3
	    a1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**5/dt
	    a2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**5/dt
	    a3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**5/dt
	    b1=(1.d0/nu2**4-1.d0/nu3**4)/(nu2*nu3)**3/dt
	    b2=(1.d0/nu3**4-1.d0/nu1**4)/(nu1*nu3)**3/dt
	    b3=(1.d0/nu1**4-1.d0/nu2**4)/(nu1*nu2)**3/dt
	    c1=(1.d0/nu3**2-1.d0/nu2**2)/(nu2*nu3)**3/dt
	    c2=(1.d0/nu1**2-1.d0/nu3**2)/(nu1*nu3)**3/dt
	    c3=(1.d0/nu2**2-1.d0/nu1**2)/(nu1*nu2)**3/dt
	    x3=0.d0
	    x5=0.d0
	    x7=0.d0
	    do i_nu=mesh(i-1)+key,mesh(i+1),key
	      v=pi2_bet*i_nu
	      sn=sin(v*tau_mesh(i_tau))
	      v=abs(v-ref)
	      v2=v*v
	      v3=v2*v
	      v5=v3*v2
	      v7=v5*v2
	      x3=x3+sn/v3
	      x5=x5+sn/v5
	      x7=x7+sn/v7
	    enddo
	    tran(i1,i_tau)=tran(i1,i_tau)+(x3*a1+x5*b1+x7*c1)*bti
	    tran(i2,i_tau)=tran(i2,i_tau)+(x3*a2+x5*b2+x7*c2)*bti
	    tran(i3,i_tau)=tran(i3,i_tau)+(x3*a3+x5*b3+x7*c3)*bti
	  enddo  !! over i
	enddo  !! over i_tau
	deallocate(mesh)
      end