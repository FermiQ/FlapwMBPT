      subroutine output_boson_dir(bos,ch)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*6, intent(in) :: ch
      real*8, intent(in) :: bos(0:n_nu,n_phead_dir)
      integer :: k,i_len,i_nu
      character*4 :: ch0
      if(ubi=='dft') ch0='.dft'
      if(ubi==' hf') ch0='.hf '
      if(ubi==' qp') ch0='.qp '
      if(ubi==' gw') ch0='.gw '
      if(ubi=='psi') ch0='.psi'
      if(ubi=='bsp') ch0='.bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_NU_DIR_'//ch//ch0)
        do i_nu=0,n_nu
          write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &                                 (bos(i_nu,k),k=1,n_phead_dir)
        enddo
        close(3)
      endif
      end
