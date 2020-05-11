      subroutine resp_lda_static(ch,n_ix0,ind_ix0)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ix0,ind_ix0(n_ix0)
      integer :: n_ch,i_len,i
      real*8 :: dq,pi2a
      real*8, allocatable :: xi0_nu(:,:),xi_nu(:,:),eps1(:),eps10(:)
      pi2a=(pi+pi)/par
      if(ch=='___00') n_ch=1
      if(ch=='___ZZ') n_ch=1
      if(ch=='00_ZZ') n_ch=3
      if(ch=='___XY') n_ch=1
      allocate(xi0_nu(n_q_long,n_ch))
      allocate(xi_nu(n_q_long,n_ch))
      allocate(eps10(n_q_long))
      allocate(eps1(n_q_long))
      xi0_nu=0.d0
      xi_nu=0.d0
      eps1=0.d0
      eps10=0.d0
c -------------------------------------------------------------------
      call chi_lda_static(ch,xi0_nu,xi_nu,n_ix0,ind_ix0,n_ch,eps1,eps10)
c -------------------------------------------------------------------
      xi_nu=xi_nu*amega
      xi0_nu=xi0_nu*amega
      if(maswrk) then
	    i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//ch//'RSP.dft')
	    do i=1,n_q_long
          dq=pi2a*step_q_long*(i-1)
	      write(3,100)dq,xi0_nu(i,1),xi_nu(i,1)
	    enddo
        close(3)
	  endif
      deallocate(xi0_nu,xi_nu,eps10,eps1)
      call timel('**** RESP_LDA_STATIC finished ******')
100	  format(3(e12.5,1x))
      end
