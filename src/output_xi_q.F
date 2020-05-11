      subroutine output_xi_q(ch,rsp)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      character*11, intent(in) :: ch
      real*8, intent(in) :: rsp(npnt)
	  integer :: i_len,iq,iq0,i
	  real*8 :: w
	  i_len=len_trim(allfile)
	  if(maswrk) then
c --------------------------------------------------------------------
	    if(ubi=='dft')open(3,file=allfile(1:i_len)//ch//'.dft')
	    if(ubi==' hf')open(3,file=allfile(1:i_len)//ch//'.hf')
	    if(ubi==' qp')open(3,file=allfile(1:i_len)//ch//'.qp')
	    if(ubi==' gw')open(3,file=allfile(1:i_len)//ch//'.gw')
	    if(ubi=='psi')open(3,file=allfile(1:i_len)//ch//'.psi')
	    if(ubi=='bsp')open(3,file=allfile(1:i_len)//ch//'.bsp')
	    do iq=1,n_line
	      iq0=kline_in_npnt(iq)
	      w=sqrt(dot_product(pnt(:,iq0),pnt(:,iq0)))
	      write(3,100)(pnt(i,iq0),i=1,3),w,rsp(iq0)
	    enddo
	    close(3)
	    write(iun,*)ch,' on all npnt k-points:'
	    write(iun,101)(rsp(i),i=1,npnt)
c --------------------------------------------------------------------	
	  endif
100	  format(4f10.4,1x,e12.5)
101	  format(5(1x,e12.5))
      end
