      subroutine v_coul_full_r(npw,v_full,ind_q,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: npw,ind_q,khbr
      real*8, intent(out) :: v_full(n_pbmt+npw,n_pbmt+npw)
      integer :: n_pbt,ind_j,j
      real*8, allocatable :: v0(:,:)
      allocate(v0(n_pbmt,n_pbmt))
      call v_coul_full_mm_r(v0,ind_q,khbr)
      v_full(1:n_pbmt,1:n_pbmt)=v0
      deallocate(v0)
      allocate(v0(n_pbmt,npw))
      call v_coul_full_mi_r(npw,v0,ind_q,khbr)
      v_full(1:n_pbmt,n_pbmt+1:n_pbmt+npw)=v0
      v_full(n_pbmt+1:n_pbmt+npw,1:n_pbmt)=transpose(v0)
      deallocate(v0)
      allocate(v0(npw,npw))
      call v_coul_full_ii_r(npw,v0,ind_q,khbr)
      v_full(n_pbmt+1:n_pbmt+npw,n_pbmt+1:n_pbmt+npw)=v0
      deallocate(v0)
      end
