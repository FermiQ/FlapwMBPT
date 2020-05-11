      subroutine output_boson_q_nu(bos,ch,nl,ck)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      use vertex_mod
      implicit none
      character*1, intent(in) :: ck
      character*6, intent(in) :: ch
      integer, intent(in) :: nl
	  complex*16, intent(in) :: bos(0:n_nu,nl)
	  real*8 :: qq,pi2a,ar,ai
	  integer :: k,i_len,i_nu,m,m1,k0
	  character*4 :: ch0
	  pi2a=(pi+pi)/par
	  if(ubi=='dft') ch0='.dft'
	  if(ubi==' hf') ch0='.hf '
	  if(ubi==' qp') ch0='.qp '
	  if(ubi==' gw') ch0='.gw '
	  if(ubi=='psi') ch0='.psi'
	  if(ubi=='bsp') ch0='.bsp'
      ar=maxval(abs(real(bos)))
      ai=maxval(abs(imag(bos)))
	  m=min(nl,5)
	  m1=min(n_nu,4)
	  if(maswrk) then
	    i_len=len_trim(allfile)
        if(ar>1.d-10) then
	      open(3,file=allfile(1:i_len)//'_NU_Q_'//ch//'_Re'//ch0)
	      do i_nu=0,n_nu
 	        write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                               (real(bos(i_nu,k)),k=1,m)
	      enddo
	      close(3)
        endif
        if(ai>1.d-10) then
	      open(3,file=allfile(1:i_len)//'_NU_Q_'//ch//'_Im'//ch0)
	      do i_nu=0,n_nu
 	        write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                               (imag(bos(i_nu,k)),k=1,m)
	      enddo
	      close(3)
        endif
        if(nl/=1.and.ar>1.d-10) then
	      open(3,file=allfile(1:i_len)//'_Q_NU_'//ch//'_Re'//ch0)
	      do k=1,nl
            if(ck=='a') then
	          k0=kline_in_npnt(k)
            else if(ck=='c') then
	          k0=kline_in_npnt_c(k)
              k0=k_a_from_c(k0)
            endif
	        qq=sqrt(dot_product(pnt(:,k0),pnt(:,k0)))*pi2a
 	        write(3,'(f11.5,5(1x,e12.5))')qq,
     &	                               (real(bos(i_nu,k)),i_nu=0,m1)
	      enddo
	      close(3)
        endif
        if(nl/=1.and.ai>1.d-10) then
	      open(3,file=allfile(1:i_len)//'_Q_NU_'//ch//'_Im'//ch0)
	      do k=1,nl
            if(ck=='a') then
	          k0=kline_in_npnt(k)
            else if(ck=='c') then
	          k0=kline_in_npnt_c(k)
              k0=k_a_from_c(k0)
            endif
	        qq=sqrt(dot_product(pnt(:,k0),pnt(:,k0)))*pi2a
 	        write(3,'(f11.5,5(1x,e12.5))')qq,
     &	                               (imag(bos(i_nu,k)),i_nu=0,m1)
	      enddo
	      close(3)
        endif
	  endif
      end
