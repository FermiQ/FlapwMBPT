      subroutine heg1_g_r_tau
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i_tau,ind_tau,it,ispin,i,ii
      real*8 :: tc,x,x2,x3
      real*8, allocatable :: coef(:,:,:),gk(:)
      g_r_tau_heg=0.d0
      if(iter_psi+iter_bsp/=0) g_r_tau_heg_c=0.d0
      do ispin=1,nspin
	  do i_tau=0,n_tau
	    call from_k_to_r_heg1(g_k_tau_heg(0,i_tau,ispin),
     &                          g_r_tau_heg(0,i_tau,ispin))
	  enddo
        do i_tau=1,3
	    call from_k_to_r_heg1(dg_k_tau_heg(0,i_tau,ispin),
     &                          dg_r_tau_heg(0,i_tau,ispin))
	  enddo
        if(iter_psi+iter_bsp/=0) then
          allocate(coef(n_tau,4,nstar_c))
          do i=1,nstar_c
            call spline_inhmg(tau_mesh,g_k_tau_heg_c(0,i,ispin),
     &                        coef(1,1,i),coef(1,2,i),coef(1,3,i),
     &                        coef(1,4,i),n_tau,1,
     &                        dg_k_tau_heg_c(1,i,ispin),
     &                        dg_k_tau_heg_c(3,i,ispin))
          enddo
          allocate(gk(nstar_c))
          do ind_tau=1,ndim3_tau
            i_tau=ndim3_tau*me_t+ind_tau-1
            do it=1,2
              if(it==2) i_tau=n_tau-i_tau
              tc=tau_mesh(i_tau)
              do i=1,n_tau
                if(tc>=tau_mesh(i-1).and.tc<=tau_mesh(i)) ii=i
              enddo
              x=tc-tau_mesh(ii-1)
              x2=x*x
              x3=x2*x
              do i=1,nstar_c
                gk(i)=coef(ii,1,i)+coef(ii,2,i)*x+coef(ii,3,i)*x2
     &               +coef(ii,4,i)*x3
              enddo
              call from_q_to_rr_heg_red(nstar_c,gk,
     &	                              g_r_tau_heg_c(1,i_tau,ispin))
            enddo
          enddo
          deallocate(coef,gk)
        endif
      enddo
      if(nproc_t/=1) then
        if(iter_psi+iter_bsp/=0) call dgop(g_r_tau_heg_c,
     &                            nrr_div_c*(n_tau+1)*nspin,'  +',
     &                            comm_t)
      endif
      end