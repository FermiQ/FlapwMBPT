      subroutine tsig_stat_tau_to_nu(t_stat,iq,npb,nkq,x,i_nu,it,trn)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: iq,npb,nkq,i_nu,it
      real*8, intent(in) :: trn(0:n_tau/2,n_nu2_max,2)
      complex*16, intent(in) :: t_stat(n_pbtot_red,nbndf_bnd,nqdiv_c,
     &                                 2,0:n_tau/2)
      complex*16, intent(inout) :: x(npb,nkq)
      integer :: i_tau,n
      complex*16 :: cc
      complex*16, allocatable :: aq(:,:,:,:),t_nu(:,:)
      allocate(aq(npb,nkq,0:n_tau/2,2))
      allocate(t_nu(npb,nkq))
      do i_tau=0,n_tau/2
        aq(:,:,i_tau,1)=t_stat(1:npb,1:nkq,iq,1,i_tau)
     &                 +t_stat(1:npb,1:nkq,iq,2,i_tau)
        aq(:,:,i_tau,2)=t_stat(1:npb,1:nkq,iq,1,i_tau)
     &                 -t_stat(1:npb,1:nkq,iq,2,i_tau)
      enddo
      n=2*npb*nkq
      call dgemm('n','n',n,1,n_tau/2+1,1.d0,aq(1,1,0,1),n,trn(0,i_nu,1),
     &           n_tau/2+1,0.d0,t_nu,n)
      x=x+t_nu
      call dgemm('n','n',n,1,n_tau/2+1,1.d0,aq(1,1,0,2),n,trn(0,i_nu,2),
     &           n_tau/2+1,0.d0,t_nu,n)
      cc=(0.d0,1.d0)
      if(it==1) cc=-cc
      x=x+cc*t_nu
      deallocate(aq,t_nu)
      end