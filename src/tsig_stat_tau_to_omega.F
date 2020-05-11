      subroutine tsig_stat_tau_to_omega(t_stat,iq,npb,nkq,x,i_omega,it,
     &                                  trn)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: iq,npb,nkq,i_omega,it
      real*8, intent(in) :: trn(0:n_tau/2,n_omega1_max,4)
      complex*16, intent(in) :: t_stat(n_pbtot_red,nbndf_bnd,nqdiv_c,
     &                                 2,0:n_tau/2)
      complex*16, intent(inout) :: x(npb,nkq)
      integer :: i_tau,n,i0
      complex*16 :: cc
      complex*16, allocatable :: aq(:,:,:,:),t_omega(:,:)
      allocate(aq(npb,nkq,0:n_tau/2,2))
      allocate(t_omega(npb,nkq))
      do i_tau=0,n_tau/2
        aq(:,:,i_tau,1)=t_stat(1:npb,1:nkq,iq,1,i_tau)
     &                 +t_stat(1:npb,1:nkq,iq,2,i_tau)
        aq(:,:,i_tau,2)=t_stat(1:npb,1:nkq,iq,1,i_tau)
     &                 -t_stat(1:npb,1:nkq,iq,2,i_tau)
      enddo
      i0=(it-1)*2+1
      n=2*npb*nkq
      call dgemm('n','n',n,1,n_tau/2+1,1.d0,aq(1,1,0,1),n,
     &           trn(0,i_omega,i0),n_tau/2+1,0.d0,t_omega,n)
      x=x+t_omega
      call dgemm('n','n',n,1,n_tau/2+1,1.d0,aq(1,1,0,2),n,
     &           trn(0,i_omega,i0+1),n_tau/2+1,0.d0,t_omega,n)
      cc=(0.d0,1.d0)
      x=x+cc*t_omega
      deallocate(aq,t_omega)
      end