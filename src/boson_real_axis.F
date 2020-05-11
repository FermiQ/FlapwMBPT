      subroutine boson_real_axis(ch,e,nl,loge,key)
c     loge = T if E or E1 - we subtract and add 1 in this case
c     key = 1 - for functions with 1/nu asymptotic
c     key = 2 - for functions with 1/nu^2 asymptotic      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      character*6, intent(in) :: ch
      logical, intent(in) :: loge
      integer, intent(in) :: nl,key
      complex*16, intent(in) :: e(0:n_nu,nl)
	character*9 :: ch2,ch3
	integer :: i_len,n,i,nk,nrax,iq
	real*8 :: omeg,ev,h_e_rax
	real*8, allocatable :: wrx(:)
	complex*16, allocatable :: s(:),sr(:,:)
	ev=evolt/2.d0
	i_len=len_trim(allfile)
c ------ Real frequency mesh --------------------------	
	nrax=nrax_chi
      allocate(wrx(0:nrax))
      h_e_rax=sqrt(freq_chi)/nrax
      do i=0,nrax
        wrx(i)=(i*h_e_rax)**2
      enddo
c ------------------------------------------------------
	allocate(sr(0:nrax,nl))
	allocate(s(0:n_nu))  
	do iq=1,nl
        s=e(:,iq)
        if(loge) s=s-1.d0
        call pade_boson_all(s,wrx,nrax,sr(0,iq),key)
        if(loge) sr(:,iq)=sr(:,iq)+1.d0
	enddo  !! over iq
	deallocate(s)
	nk=min(5,nl)
	ch2=ch//'_Re_'
	ch3=ch//'_Im_'
	if(maswrk) then
	  if(ubi=='dft') open(3,file=allfile(1:i_len)//ch2//'dft.rax')
	  if(ubi==' hf') open(3,file=allfile(1:i_len)//ch2//'hf.rax')
	  if(ubi==' qp') open(3,file=allfile(1:i_len)//ch2//'qp.rax')
	  if(ubi==' gw') open(3,file=allfile(1:i_len)//ch2//'gw.rax')
	  if(ubi=='psi') open(3,file=allfile(1:i_len)//ch2//'psi.rax')
	  if(ubi=='bsp') open(3,file=allfile(1:i_len)//ch2//'bsp.rax')
	  do i=0,nrax
	    omeg=wrx(i)
	    write(3,100)omeg*ev,(real(sr(i,n)),n=1,nk)
        enddo
        close(3)
	  if(ubi=='dft') open(3,file=allfile(1:i_len)//ch3//'dft.rax')
	  if(ubi==' hf') open(3,file=allfile(1:i_len)//ch3//'hf.rax')
	  if(ubi==' qp') open(3,file=allfile(1:i_len)//ch3//'qp.rax')
	  if(ubi==' gw') open(3,file=allfile(1:i_len)//ch3//'gw.rax')
	  if(ubi=='psi') open(3,file=allfile(1:i_len)//ch3//'psi.rax')
	  if(ubi=='bsp') open(3,file=allfile(1:i_len)//ch3//'bsp.rax')
	  do i=0,nrax
	    omeg=wrx(i)
	    write(3,100)omeg*ev,(imag(sr(i,n)),n=1,nk)
        enddo
	  close(3)
	endif
	deallocate(sr,wrx)
100	format(f11.5,1x,5(e12.5,1x))
      end
