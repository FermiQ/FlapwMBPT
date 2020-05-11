      subroutine w_tau_prepare_red(w_tau)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  real*8, intent(out) :: w_tau(n_pbtot_red,n_pbtot_red,
     &                             ndim3_tau,nqdiv_c)
	  integer :: k,ig,ind_q,iq,n,kst,ind_nu,i_nu,nn,n0,ka
	  real*8, allocatable :: wnu(:,:,:)
	  complex*16, allocatable :: tmp(:,:)
      n=n_pbtot_red
      nn=n*n
      allocate(wnu(n,n,0:n_nu))
c ----------------------------------------------------------------
      allocate(tmp(n,n))
      w_tau=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        k=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(k)
        n0=n_pbmt_red+nplwgw_red(ka)
	    do kst=1,k_star_c(k)
	      iq=k_list_c(kst,k)
	      ig=k_group_c(iq)
	      wnu=0.d0
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            tmp=w_red_q(:,:,ind_nu,ka)
	        call sym_w_red(iq,k,tmp)
	        call pack_hermit(tmp,wnu(1,1,i_nu),n,n,n,0.d0,1.d0)
          enddo  !! over ind_nu
	      if(nproc_t/=1) call dgop(wnu,nn*(n_nu+1),'  +',comm_t)
          call nua_to_tau_ba(wnu,n,w_tau(1,1,1,iq),n,n0)
        enddo   !! over kst
      enddo   !! over ind_q
      deallocate(tmp,wnu)
	  if(nproc_k/=1) call dgop(w_tau,nn*nqdiv_c*ndim3_tau,'  +',
     &                         comm_k)
      end
