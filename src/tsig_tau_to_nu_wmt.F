      subroutine tsig_tau_to_nu_wmt(t_pw,ind_omega,nomnu,ispin,numtn,
     &                              np,n1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,nomnu,ispin,numtn,np,n1
      complex*16, intent(inout) :: t_pw(np,n1,numtn,2)
      integer :: it,i_tau,j,i,i1_tau
      complex*16 :: cc
      complex*16, allocatable :: vxt(:,:,:,:),aq(:,:,:,:),t_nu(:,:,:)
c -------- Transform Vrt(tau;omega) ---- > Vrt(nu;omega) ------------
      allocate(vxt(np,n1,0:n_tau,2))
      do it=1,2
        do i_tau=0,n_tau
          vxt(:,:,i_tau,it)=t_pw(:,:,i_tau+1,it)
        enddo
      enddo
      allocate(aq(np,n1,0:n_tau/2,4))
      allocate(t_nu(np,n1,nomnu))
      do it=1,2
        cc=(-1)**it*(0.d0,1.d0)
        aq=(0.d0,0.d0)
        do i_tau=0,n_tau/2
          i1_tau=n_tau-i_tau
          do j=1,n1
            do i=1,np
              aq(i,j,i_tau,1)=vxt(i,j,i_tau,3-it)+vxt(i,j,i1_tau,3-it)
              aq(i,j,i_tau,2)=vxt(i,j,i_tau,3-it)-vxt(i,j,i1_tau,3-it)
              aq(i,j,i_tau,2)=cc*aq(i,j,i_tau,2)
              aq(i,j,i_tau,3)=vxt(i,j,i_tau,it)+vxt(i,j,i1_tau,it)
              aq(i,j,i_tau,4)=vxt(i,j,i_tau,it)-vxt(i,j,i1_tau,it)
              aq(i,j,i_tau,4)=cc*aq(i,j,i_tau,4)
            enddo
          enddo
        enddo
        call from_tau_om_to_nu_om(2*np*n1,ind_omega,nomnu,t_nu,aq)
        do i=1,nomnu
          t_pw(:,:,i,it)=t_nu(:,:,i)
        enddo
      enddo
      deallocate(vxt,aq,t_nu)
      end
