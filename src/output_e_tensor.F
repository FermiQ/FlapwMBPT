      subroutine output_e_tensor(e,ch)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      complex*16, intent(in) :: e(3,3,0:n_nu)
      integer :: i_len,i_nu
      character*4 :: ch0
      if(ubi=='dft') ch0='.dft'
      if(ubi==' hf') ch0='.hf '
      if(ubi==' qp') ch0='.qp '
      if(ubi==' gw') ch0='.gw '
      if(ubi=='psi') ch0='.psi'
      if(ubi=='bsp') ch0='.bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_NU_Q_'//ch//'_Tensor_Re'//ch0)
        do i_nu=0,n_nu
          write(3,'(f11.5,6(1x,f11.6))')w_nu(i_nu)*evolt/2,
     &                                 real(e(1,1,i_nu)),
     &                                 real(e(2,2,i_nu)),
     &                                 real(e(3,3,i_nu)),
     &                                 real(e(1,2,i_nu)),
     &                                 real(e(1,3,i_nu)),
     &                                 real(e(2,3,i_nu))
        enddo
        close(3)
      endif
      end
