      subroutine output_p_tensor(p)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      complex*16, intent(in) :: p(3,3,ndim3_nu)
      integer :: i_len,i_nu,ind_nu
      character*4 :: ch0
      real*8, allocatable :: p0(:,:,:)
      allocate(p0(3,3,0:n_nu))
      p0=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        p0(:,:,i_nu)=p(:,:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(p0,9*(n_nu+1),'  +',comm_t)
      if(ubi=='dft') ch0='.dft'
      if(ubi==' hf') ch0='.hf '
      if(ubi==' qp') ch0='.qp '
      if(ubi==' gw') ch0='.gw '
      if(ubi=='psi') ch0='.psi'
      if(ubi=='bsp') ch0='.bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_NU_Q_P'//'_Tensor_Re'//ch0)
        do i_nu=0,n_nu
          write(3,'(f11.5,6(1x,f11.6))')w_nu(i_nu)*evolt/2,
     &                                 p0(1,1,i_nu),p0(2,2,i_nu),
     &                                 p0(3,3,i_nu),p0(1,2,i_nu),
     &                                 p0(1,3,i_nu),p0(2,3,i_nu)
        enddo
        close(3)
      endif
      deallocate(p0)
      end


      subroutine output_p_tensor_r(p)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: p(3,3,ndim3_nu)
      integer :: i_len,i_nu,ind_nu
      character*4 :: ch0
      real*8, allocatable :: p0(:,:,:)
      allocate(p0(3,3,0:n_nu))
      p0=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        p0(:,:,i_nu)=p(:,:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(p0,9*(n_nu+1),'  +',comm_t)
      if(ubi=='dft') ch0='.dft'
      if(ubi==' hf') ch0='.hf '
      if(ubi==' qp') ch0='.qp '
      if(ubi==' gw') ch0='.gw '
      if(ubi=='psi') ch0='.psi'
      if(ubi=='bsp') ch0='.bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_NU_Q_P'//'_Tensor_Re'//ch0)
        do i_nu=0,n_nu
          write(3,'(f11.5,6(1x,f11.6))')w_nu(i_nu)*evolt/2,
     &                                 p0(1,1,i_nu),p0(2,2,i_nu),
     &                                 p0(3,3,i_nu),p0(1,2,i_nu),
     &                                 p0(1,3,i_nu),p0(2,3,i_nu)
        enddo
        close(3)
      endif
      deallocate(p0)
      end
