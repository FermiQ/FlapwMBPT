      subroutine sym_w_lim(iq,iq0,w_0)
      use atom_mod
      use solid_mod
      implicit none
      integer, intent(in) :: iq,iq0
      complex*16, intent(inout) :: w_0(n_pb_lim,n_pb_lim)
      integer :: j,ig
      complex*16, allocatable :: tmp(:)
c -------- Copy+Symmetrization of W into another k-point ---------------
      ig=k_group(iq)
      if(iq/=iq0) then
        call sym_w_mt_lim(ig,pnt(1,iq),w_0,n_pb_lim)
        allocate(tmp(n_pb_lim-1))
        tmp=w_0(1:n_pb_lim-1,n_pb_lim)
        call sym_w_mi_lim(iq,ig,tmp)
        w_0(1:n_pb_lim-1,n_pb_lim)=tmp
        deallocate(tmp)
      endif
c ------- Obtaining  W_im from W_mi -------------------------------
      do j=1,n_pb_lim-1
        w_0(n_pb_lim,j)=conjg(w_0(j,n_pb_lim))
      enddo
      end
