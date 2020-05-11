      subroutine output_xi_cmpl(ch,f,cnu,nfreq,key)
c     key = 0 - for Xi 
c     key = 1 - for Ixc
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*9, intent(in) :: ch
      integer, intent(in) :: nfreq,key
      complex*16, intent(in) :: cnu(nfreq)
	complex*16, intent(in) :: f(nfreq,n_line)
	integer :: i_len,n,i,nk
	real*8 :: ev,ev1
	ev=evolt/2.d0
	ev1=ev
	if(key==0) ev1=1.d0/ev
	i_len=len_trim(allfile)
	nk=min(4,n_line)
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//ch//'_Re')
	  write(3,*)' First two columns - Re and Im parts of Frequency'
	  do i=1,nfreq
	    write(3,100)cnu(i)*ev,(real(f(i,n)*ev1),n=1,nk)
        enddo
        close(3)
	  open(3,file=allfile(1:i_len)//ch//'_Im')
	  write(3,*)' First two columns - Re and Im parts of Frequency'
	  do i=1,nfreq
	    write(3,100)cnu(i)*ev,(imag(f(i,n)*ev1),n=1,nk)
        enddo
        close(3)
	endif
100	format(2(f11.5,1x),4(e12.5,1x))
      end