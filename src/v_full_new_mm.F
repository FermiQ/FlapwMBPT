      subroutine v_full_new_mm(v_full,k,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,khbr
      complex*16, intent(out) :: v_full(n_pbmt,n_pbmt)
      integer :: j,i
      real*8 :: q2a
      complex*16, allocatable :: ss(:,:)
      allocate(ss(nlbs,nlbs))
      call str_con(ss,pnt(1,k),khbr,nlbs,nlbs,'s','s')
      call v_coul_new_mt(v_full,ss,khbr)
      deallocate(ss)
c ---------- Correction terms --------------------------------------
      if(khbr==0.and.k==1) call v_mm_nonsing(v_full)
      end
