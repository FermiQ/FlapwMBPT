      subroutine from_nu_to_rax_mt(ch,f)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      complex*16, intent(inout) :: f(0:n_nu,nsort)
      integer :: i,nrax,isort
      real*8 :: h_xi_rax
      real*8, allocatable :: wrx(:)
      complex*16, allocatable :: xi_rax(:,:,:)
      nrax=nrax_chi
      allocate(xi_rax(0:nrax,nsort,2))
      xi_rax=(0.d0,0.d0)
      allocate(wrx(0:nrax))
      h_xi_rax=sqrt(freq_chi)/nrax
      do i=0,nrax
        wrx(i)=(i*h_xi_rax)**2
      enddo
      do isort=1,nsort
        call pade_boson_all(f(0,isort),wrx,nrax,xi_rax(0,isort,1),2)
      enddo  !! over isort
      if(ch=='+-') then
        f=conjg(f)
        do isort=1,nsort
          call pade_boson_all(f(0,isort),wrx,nrax,xi_rax(0,isort,2),2)
        enddo  !! over isort
        f=conjg(f)
      endif
      call output_xi_rax_mt(ch,xi_rax,wrx,nrax)
      deallocate(xi_rax,wrx)
      end
