      subroutine optics
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      character*1 :: mode
      integer :: nrax,i,inu
      real*8 :: h_e_rax,f
      complex*16 :: tc(3,3)
      real*8, allocatable :: sig_r(:,:,:),eloss(:,:,:)
      complex*16, allocatable :: wrx(:),e(:,:,:)
      if(irel/=0) return
      if(ubi=='dft') mode=opt_mode(1:1)
      if(ubi==' hf') mode=opt_mode(2:2)
      if(ubi==' gw') mode=opt_mode(3:3)
      if(ubi==' qp') mode=opt_mode(4:4)
      if(ubi=='psi') mode=opt_mode(5:5)
      if(ubi=='bsp') mode=opt_mode(6:6)
      if(mode=='0') return
c ------ Real frequency mesh --------------------------	
	  nrax=nrax_chi
      allocate(wrx(0:nrax))
      h_e_rax=sqrt(freq_chi)/nrax
      do i=0,nrax
        f=(i*h_e_rax)**2
        wrx(i)=dcmplx(f,e_small_bos)
      enddo
c ------------------------------------------------------
      allocate(e(3,3,0:nrax))
      if(ubi=='dft') then
c        if(.not.opt_loc_fld) call optics_hf_dir(nrax,wrx,e,mode)
c        if(.not.opt_loc_fld) call optics_no_loc_fld(nrax,wrx,e)
c        if(opt_loc_fld) call optics_loc_fld(nrax,wrx,e,mode)
c        call optics_gen_no_lcf(nrax,wrx,e)
        call optics_dir_no_lcf(nrax,wrx,e)
      else if(ubi==' hf') then
        if(.not.opt_loc_fld) call optics_hf_dir(nrax,wrx,e,mode)
c        call optics_hf(nrax,wrx,e)
      else
        call optics_gen_no_lcf(nrax,wrx,e)
      endif
      allocate(sig_r(3,3,0:nrax))
      do inu=0,nrax
        sig_r(:,:,inu)=real(wrx(inu))*imag(e(:,:,inu))/4.d0/pi
      enddo
      allocate(eloss(3,3,0:nrax))
      do inu=0,nrax
        tc=e(:,:,inu)
        call invers_z(3,tc,3)
        eloss(:,:,inu)=-imag(tc)
      enddo
      call output_optics(e,sig_r,eloss,nrax,wrx)
      deallocate(wrx,e,sig_r,eloss)
      call timel('***** OPTICS finished **************')
      end
