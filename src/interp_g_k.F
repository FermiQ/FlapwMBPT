      subroutine interp_g_k(f,x,y)
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: f(0:n_k_heg),x
	real*8, intent(out) :: y
	integer :: i0
	real*8 :: xm(0:1),t
	if(x>k_infinum) then
	  y=0.d0
	  return
	endif
	t=x_fermi+asinh((x-k_fermi)/alf_k_heg)
	i0=t/step_k_x
	xm(0)=i0*step_k_x
	xm(1)=xm(0)+step_k_x
      call fint(f(i0),xm,2,t,y)
      end
      
      subroutine interp_g_k_c(f,x,y)
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: x
	complex*16, intent(in) :: f(0:n_k_heg)
	complex*16, intent(out) :: y
	integer :: i0
	real*8 :: xm(0:1),t,y1,y2
	if(x>k_infinum) then
	  y=(0.d0,0.d0)
	  return
	endif
	t=x_fermi+asinh((x-k_fermi)/alf_k_heg)
	i0=t/step_k_x
	if(i0==n_k_heg) then
	  y=f(n_k_heg)
	else
	  xm(0)=i0*step_k_x
	  xm(1)=xm(0)+step_k_x
        call fint(real(f(i0:i0+1)),xm,2,t,y1)
        call fint(imag(f(i0:i0+1)),xm,2,t,y2)
        y=dcmplx(y1,y2)
      endif
      end
      
      subroutine interp_k_x_c(f,k,x,y)
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: k,x
	complex*16, intent(in) :: f(0:n_k_heg,0:n_k_heg)
	complex*16, intent(out) :: y
	integer :: ik,ix,ik1,ix1
	real*8 :: t,dk,dx,s,s0,s1,s2
	if(k>k_infinum) then
	  y=(0.d0,0.d0)
	  return
	endif
	t=x_fermi+asinh((k-k_fermi)/alf_k_heg)
	ik=idnint(t/step_k_x)  ! closest k-point
	dk=t-ik*step_k_x
	if(dk>1.d-14) then
	  ik1=ik+1
	else if(dk<-1.d-14) then
	  ik1=ik-1
	else
	  ik1=ik
	endif
	t=x+1.d0
	ix=idnint(t/step_x_x)  ! closest x-point
	dx=t-ix*step_x_x
	if(dx>1.d-12) then
	  ix1=ix+1
	else if(dx<-1.d-12) then
	  ix1=ix-1
	else
	  ix1=ix
	endif
	s=step_k_x*step_x_x
	s1=step_x_x*abs(dk)
	s2=step_k_x*abs(dx)
	s0=s-s1-s2
      y=(s0*f(ik,ix)+s1*f(ik1,ix)+s2*f(ik,ix1))/s
      end
      
      subroutine interp_triangle_index(k,x,i1,i2,i3,c1,c2,c3)
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	real*8, intent(in) :: k,x
	integer, intent(out) :: i1,i2,i3
	real*8, intent(out) :: c1,c2,c3
	integer :: ik,ix,ik1,ix1
	real*8 :: t,dk,dx,s,s0,s1,s2
	if(k>k_infinum) then
	  ik=n_k_heg
	  ik1=ik
	  dk=0.d0
	else
	  t=x_fermi+asinh((k-k_fermi)/alf_k_heg)
	  ik=idnint(t/step_k_x)  ! closest k-point
	  dk=t-ik*step_k_x
	  if(dk>1.d-14) then
	    ik1=ik+1
	  else if(dk<-1.d-14) then
	    ik1=ik-1
	  else
	    ik1=ik
	  endif
	endif
	t=x+1.d0
	ix=idnint(t/step_x_x)  ! closest x-point
	dx=t-ix*step_x_x
	if(dx>1.d-12) then
	  ix1=ix+1
	else if(dx<-1.d-12) then
	  ix1=ix-1
	else
	  ix1=ix
	endif
	s=step_k_x*step_x_x
	s1=step_x_x*abs(dk)
	s2=step_k_x*abs(dx)
	s0=s-s1-s2
	i1=(n_k_heg+1)*ix+ik+1
	i2=(n_k_heg+1)*ix+ik1+1
	i3=(n_k_heg+1)*ix1+ik+1
	c1=s0/s
	c2=s1/s
	c3=s2/s
      end