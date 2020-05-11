      subroutine sym_p_red(iq,ig,w_0)
      use atom_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: iq,ig
      complex*16, intent(inout) :: w_0(n_pbtot_red,n_pbtot_red)
      integer :: j,i,iq0
      complex*16, allocatable :: tmp(:,:)
c -------- Copy+Symmetrization of W into another k-point ---------------
      iq0=ip_k(iq,ig)
      call sym_w_mt_red(ig,pnt(1,iq),w_0,n_pbtot_red)
      allocate(tmp(n_pbmt_red,nplw_gw_red))
      tmp=w_0(1:n_pbmt_red,n_pbmt_red+1:n_pbtot_red)
      call sym_p_mi_red(iq,ig,tmp)
      w_0(1:n_pbmt_red,n_pbmt_red+1:n_pbtot_red)=tmp
      deallocate(tmp)
      allocate(tmp(nplw_gw_red,nplw_gw_red))
      tmp=w_0(n_pbmt_red+1:n_pbtot_red,n_pbmt_red+1:n_pbtot_red)
      call sym_p_ii_red(iq,ig,tmp,nplwgw_red(iq0),nplw_gw_red)
      w_0(n_pbmt_red+1:n_pbtot_red,n_pbmt_red+1:n_pbtot_red)=tmp
      deallocate(tmp)
c ------- Obtaining  W_im from W_mi -------------------------------
      do j=1,n_pbmt_red
        do i=n_pbmt_red+1,n_pbtot_red
          w_0(i,j)=conjg(w_0(j,i))
        enddo
      enddo
      end
