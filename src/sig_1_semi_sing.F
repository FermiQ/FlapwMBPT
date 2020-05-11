      subroutine sig_1_semi_sing(t_sing,g_tau,s2,sig_semi,nb)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: s2,nb
      real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau)
      complex*16, intent(in) :: t_sing(nbndf_bnd,2,ndim3_tau)
      complex*16, intent(inout) :: sig_semi(nbndf_bnd,nbndf_bnd,2,
     &                                      ndim3_tau)
      integer :: it,ind_nu,i_nu,ind_tau,i,j
      real*8 :: w
      complex*16, allocatable :: t_nu(:,:,:),wt(:,:,:),wtt(:,:,:),
     &                           gg(:,:),tt(:)
      allocate(t_nu(nbndf_bnd,ndim3_nu,2))
      call tau_to_nu_full(t_nu,t_sing,nbndf_bnd)
      allocate(wt(nb,ndim3_nu,2))
      do it=1,2
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          w=ws_head_nu(i_nu)*q2aver_c
          wt(:,ind_nu,it)=w*t_nu(1:nb,ind_nu,it)
        enddo
      enddo
      deallocate(t_nu)
      allocate(wtt(nb,2,ndim3_tau))
      call nu_to_tau_cc(wt,wtt,nb)
      deallocate(wt)
      allocate(gg(nb,nb))
      allocate(tt(nb))
      w=1.d0/sqrt(amega)/nqdiv_c
      do ind_tau=1,ndim3_tau
        do it=1,2
          call ferm_unpack_tau(gg,g_tau(1,1,1,ind_tau),nb,nb,
     &                         nbndf_bnd,it)
          tt=(0.d0,0.d0)
          do j=1,nb
            do i=1,nb
              tt(i)=tt(i)+gg(i,j)*wtt(j,it,ind_tau)
            enddo
          enddo
          do i=1,nb
            sig_semi(i,s2,it,ind_tau)=sig_semi(i,s2,it,ind_tau)+w*tt(i)
          enddo
        enddo
      enddo
      deallocate(wtt,gg,tt)
      end
