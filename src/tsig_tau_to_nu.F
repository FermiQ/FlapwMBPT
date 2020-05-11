      subroutine tsig_tau_to_nu(t_pw,ind_omega,k1,nomnu,ispin,numtn)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_omega,k1,nomnu,ispin,numtn
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c,
     &                                  numtn,2)
      integer :: nk,kq,kq0,nkq,n2,it,i_tau,kqa0,iq,iq0,iqa0,nk0,j,
     &           i,i1_tau
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: vxt(:,:,:,:),aq(:,:,:,:),t_nu(:,:,:)
      n2=nbndf_bnd**2
c -------- Transform Vrt(tau;omega) ---- > Vrt(nu;omega) ------------
	  do iq=1,nqdiv_c
	    iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
	    nk=nplwgw_red(iqa0)
        nk0=n_pbmt_red+nk
	    v=pnt_c(:,k1)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K'-Q
	    kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kqa0,ispin)
	    allocate(vxt(nk0,nkq,0:n_tau,2))
        do it=1,2
	      do i_tau=0,n_tau
            vxt(:,:,i_tau,it)=t_pw(1:nk0,1:nkq,iq,i_tau+1,it)
          enddo
        enddo
        allocate(aq(nk0,nkq,0:n_tau/2,4))
        allocate(t_nu(nk0,nkq,nomnu))
        do it=1,2
          cc=(-1)**it*(0.d0,1.d0)
          aq=(0.d0,0.d0)
          do i_tau=0,n_tau/2
            i1_tau=n_tau-i_tau
            do j=1,nkq
              do i=1,nk0
                aq(i,j,i_tau,1)=vxt(i,j,i_tau,3-it)
     &                         +vxt(i,j,i1_tau,3-it)
                aq(i,j,i_tau,2)=vxt(i,j,i_tau,3-it)
     &                         -vxt(i,j,i1_tau,3-it)
                aq(i,j,i_tau,2)=cc*aq(i,j,i_tau,2)
                aq(i,j,i_tau,3)=vxt(i,j,i_tau,it)
     &                         +vxt(i,j,i1_tau,it)
                aq(i,j,i_tau,4)=vxt(i,j,i_tau,it)
     &                         -vxt(i,j,i1_tau,it)
                aq(i,j,i_tau,4)=cc*aq(i,j,i_tau,4)
              enddo
            enddo
          enddo
          call from_tau_om_to_nu_om(2*nk0*nkq,ind_omega,nomnu,t_nu,aq)
          do i=1,nomnu
            t_pw(1:nk0,1:nkq,iq,i,it)=t_nu(:,:,i)
          enddo
        enddo
        deallocate(vxt,aq,t_nu)
      enddo   !! over iq
      end
