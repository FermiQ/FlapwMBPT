      subroutine tsig0_tau_to_nu(t_sing,ind_omega,k1,nomnu,ispin,numtn)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,nomnu,ispin,numtn
      complex*16, intent(inout) :: t_sing(nbndf_bnd,numtn,2)
      integer :: n,ka,it,i_tau,i1_tau,i
      complex*16 :: cc
      complex*16, allocatable :: vxt(:,:,:),aq(:,:,:),t_nu(:,:)
      ka=k_a_from_c(k1)
      n=n_low_bnd(ka,ispin)
      allocate(vxt(n,0:n_tau,2))
      do it=1,2
        do i_tau=0,n_tau
          vxt(:,i_tau,it)=t_sing(1:n,i_tau+1,it)
        enddo
      enddo
      allocate(aq(n,0:n_tau/2,4))
      allocate(t_nu(n,nomnu))
      do it=1,2
        cc=(-1)**it*(0.d0,1.d0)
        aq=(0.d0,0.d0)
        do i_tau=0,n_tau/2
          i1_tau=n_tau-i_tau
          do i=1,n
            aq(i,i_tau,1)=vxt(i,i_tau,3-it)+vxt(i,i1_tau,3-it)
            aq(i,i_tau,2)=vxt(i,i_tau,3-it)-vxt(i,i1_tau,3-it)
            aq(i,i_tau,2)=cc*aq(i,i_tau,2)
            aq(i,i_tau,3)=vxt(i,i_tau,it)+vxt(i,i1_tau,it)
            aq(i,i_tau,4)=vxt(i,i_tau,it)-vxt(i,i1_tau,it)
            aq(i,i_tau,4)=cc*aq(i,i_tau,4)
          enddo
        enddo
        call from_tau_om_to_nu_om(2*n,ind_omega,nomnu,t_nu,aq)
        do i=1,nomnu
          t_sing(1:n,i,it)=t_nu(:,i)
        enddo
      enddo
      deallocate(vxt,aq,t_nu)
      end
