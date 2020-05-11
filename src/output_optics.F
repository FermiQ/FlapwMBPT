      subroutine output_optics(e,sig_r,eloss,nrax,wrx)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
      integer, intent(in) :: nrax
      real*8, intent(in) :: sig_r(3,3,0:nrax),eloss(3,3,0:nrax)
      complex*16, intent(in) :: wrx(0:nrax),e(3,3,0:nrax)
	  integer :: i,inu,i_len
	  real*8 :: x(6)
      real*8, allocatable :: rnu(:)
	  if(.not.maswrk) return
      allocate(rnu(0:nrax))
      rnu=real(wrx)*evolt/2
	  i_len=len_trim(allfile)
	  open(3,file=allfile(1:i_len)//'_E_Optics_Re')
      do inu=0,nrax
        x(1)=real(e(1,1,inu))
        x(2)=real(e(2,2,inu))
        x(3)=real(e(3,3,inu))
        x(4)=real(e(1,2,inu))
        x(5)=real(e(1,3,inu))
        x(6)=real(e(2,3,inu))
        write(3,'(f8.5,6(1x,f10.6))')rnu(inu),(x(i),i=1,6)
      enddo
      close(3)
	  open(3,file=allfile(1:i_len)//'_E_Optics_Im')
      do inu=0,nrax
        x(1)=imag(e(1,1,inu))
        x(2)=imag(e(2,2,inu))
        x(3)=imag(e(3,3,inu))
        x(4)=imag(e(1,2,inu))
        x(5)=imag(e(1,3,inu))
        x(6)=imag(e(2,3,inu))
        write(3,'(f8.5,6(1x,f10.6))')rnu(inu),(x(i),i=1,6)
      enddo
      close(3)
	  open(3,file=allfile(1:i_len)//'_Cond_Optics_Re')
      do inu=0,nrax
        x(1)=sig_r(1,1,inu)
        x(2)=sig_r(2,2,inu)
        x(3)=sig_r(3,3,inu)
        x(4)=sig_r(1,2,inu)
        x(5)=sig_r(1,3,inu)
        x(6)=sig_r(2,3,inu)
        write(3,'(f8.5,6(1x,f10.7))')rnu(inu),(x(i),i=1,6)
      enddo
      close(3)
	  open(3,file=allfile(1:i_len)//'_ELoss_Optics_Re')
      do inu=0,nrax
        x(1)=eloss(1,1,inu)
        x(2)=eloss(2,2,inu)
        x(3)=eloss(3,3,inu)
        x(4)=eloss(1,2,inu)
        x(5)=eloss(1,3,inu)
        x(6)=eloss(2,3,inu)
        write(3,'(f8.5,6(1x,f10.7))')rnu(inu),(x(i),i=1,6)
      enddo
      close(3)
      deallocate(rnu)
      end
