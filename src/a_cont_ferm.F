      subroutine a_cont_ferm(nrax,erax,f_in,f_out,esm)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nrax
      real*8, intent(in) :: erax(nrax),esm
      complex*16, intent(in) :: f_in(0:n_omega)
      complex*16, intent(out) :: f_out(nrax)
      integer :: i
      real*8 :: chk
      complex*16 :: om0
      complex*16, allocatable :: a(:)
      chk=maxval(abs(f_in))
      if(chk<1.d-12) then
        f_out=(0.d0,0.d0)
      else
        allocate(a(0:n_omega))
        call pade_fermion_all(f_in,1,1,a)
        do i=1,nrax
          om0=dcmplx(erax(i),esm)
          call pade_fermion_all_out(1,1,a,om0,f_out(i))
        enddo
        deallocate(a)
      endif
      end
