      subroutine resp0_lda_q1_pw(ch)
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
      integer :: i_len,i
      real*8 :: h_xi_rax,ev
      real*8, allocatable :: wrx(:)
      complex*16, allocatable :: xi0_nu(:)
      ev=evolt/2.d0
      allocate(xi0_nu(0:nrax_chi))
      allocate(wrx(0:nrax_chi))
      h_xi_rax=sqrt(freq_chi)/nrax_chi
      do i=0,nrax_chi
        wrx(i)=(i*h_xi_rax)**2
      enddo
c ------------------------------------------------------------------
      call chi0_lda_dynamic_pw(ch,wrx,xi0_nu)
c ------------------------------------------------------------------
      xi0_nu=xi0_nu*amega
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//ch//'RSP0.dft')
        do i=0,nrax_chi
          write(3,100)wrx(i)*ev,xi0_nu(i)/ev
        enddo
        close(3)
      endif
      deallocate(xi0_nu,wrx)
      call timel('**** RESP0_LDA_Q1_PW finished ******')
100   format(3(e12.5,1x))
      end
