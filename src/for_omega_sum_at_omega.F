      subroutine for_omega_sum_at_omega(n_ex,n_as,tran,ncheb,ic,w0,
     &                                  w1,w2,w3,w4,w5,wg,wg1)
c     Creates transformation matrices for SUM OMEGA'_omega ---> omega
c       OMEGA'-mesh is adaptive for every omega ---------------
c       They correspond to the transform :
c                A(w) = 1/b Sum_w'  A(w';w)
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
	integer, intent(in) :: n_ex(2),n_as(3),ncheb(3),ic(4)
      real*8, intent(in) :: w0,w1,w2,w3,w4,w5,
     &                      wg(-1:n_omega_asy,n_omega_asy/2,2),
     &                      wg1(-1:n_omega_asy+2,n_omega_asy/2)
      real*8, intent(out) :: tran(n_omega1_max)
	integer :: i_omega,n,nn,i,ii
	real*8 :: bt,pi2_bet,pi_bet,v,x,ww0,wwm,a,b,c,vm,ww
      real*8, allocatable :: cheb(:),tmp(:),tr(:)
c ----------------------------------------------------------- 
	bt=1.d0/betta_t
	pi_bet=pi/betta_t
	pi2_bet=2.d0*pi/betta_t
c ---------------------- Exact Part ------------------------------	  
	do i_omega=n_as(1)+1,n_as(1)+n_ex(1)
	  tran(i_omega)=bt
	enddo
	if(n_as(2)/=0) then 
	  do i_omega=n_as(1)+n_ex(1)+n_as(2)+1,
     &	           n_as(1)+n_ex(1)+n_as(2)+n_ex(2)
	    tran(i_omega)=bt
	  enddo
	endif
	allocate(cheb(0:n_omega_asy/2-1))
	allocate(tmp(0:n_omega_asy/2-1))
c ----------------- Chebyshev's Part (Left) -----------------------
      allocate(tr(-1:n_omega_asy))
	n=ncheb(1)
	nn=n_as(1)-ic(1)+1
	tmp=0.d0
	v=w0
	ww0=sqrt(-w1)
      wwm=sqrt(-w0)
      a=2.d0*ww0*wwm
      b=ww0+wwm
      c=wwm-ww0
      vm=w1-pi2_bet
1     ww=sqrt(-v)
      x=(b-a/ww)/c
	call dcheb_recurr(x,cheb,n-1)
	tmp(0:n-1)=tmp(0:n-1)+bt*cheb(0:n-1)
	v=v+pi2_bet
	if(v-vm<1.d-12) goto 1
	call dgemm('n','n',nn,1,n,1.d0,wg(ic(1),1,1),n_omega_asy+2,tmp(0),
     &	         n_omega_asy/2,0.d0,tr(ic(1)),n_omega_asy+2)
      do i=ic(1),n_as(1)
        ii=n_as(1)-i+1
        tran(ii)=tran(ii)+tr(i)
      enddo
      deallocate(tr)
c ----------------- Chebyshev's Part (Internal) -----------------------
      if(n_as(2)/=0) then
	  n=ncheb(2)
	  nn=n_as(2)+ic(3)-ic(2)+1
	  tmp=0.d0
	  v=w2+pi2_bet
	  a=0.5d0*(w2+w3)
	  b=0.5d0*(w3-w2)
        vm=w3-pi2_bet
2	  x=(v-a)/b
	  call dcheb_recurr(x,cheb,n-1)
	  tmp(0:n-1)=tmp(0:n-1)+bt*cheb(0:n-1)
	  v=v+pi2_bet
	  if(v-vm<1.d-12) goto 2
	  call dgemm('n','n',nn,1,n,1.d0,wg1(ic(2),1),n_omega_asy+4,
     &             tmp(0),n_omega_asy/2,1.d0,
     &             tran(n_as(1)+n_ex(1)+ic(2)),n_omega1_max)
	endif
c ----------------- Chebyshev's Part (Right) -----------------------
	n=ncheb(3)
	nn=n_as(3)-ic(4)+1
	tmp=0.d0
	v=w4+pi2_bet
	ww0=sqrt(w4)
      wwm=sqrt(w5)
      a=2.d0*ww0*wwm
      b=ww0+wwm
      c=wwm-ww0
3     ww=sqrt(v)
      x=(b-a/ww)/c
	call dcheb_recurr(x,cheb,n-1)
	tmp(0:n-1)=tmp(0:n-1)+bt*cheb(0:n-1)
	v=v+pi2_bet
	if(v-w5<1.d-12) goto 3
	call dgemm('n','n',nn,1,n,1.d0,wg(ic(4),1,2),n_omega_asy+2,tmp(0),
     &	         n_omega_asy/2,1.d0,
     &           tran(n_as(1)+n_ex(1)+n_as(2)+n_ex(2)+ic(4)),
     &           n_omega_b_max)
	deallocate(cheb,tmp)
      end