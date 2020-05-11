      subroutine screened_w
      use atom_mod
      use manager_mod
	  use models_mod
	  use scalapack_mod
	  use solid_mod
	  use parallel_mod
	  use vertex_mod
	  implicit none
	  integer :: ind_k,k,n_pbt
	  complex*16, allocatable :: tmm(:,:),tmi(:,:),tii(:,:)
	  if(itest/=0.and..not.lsclp) then
	    allocate(tmm(n_pbmt,n_pbmt))
	    allocate(tmi(n_pbmt,nplw_gw))
	    allocate(tii(nplw_gw,nplw_gw))
	    if(maswrk) write(iun,*)' Symmetry check for P:'
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n_pbt=n_pbmt+nplwgw(k)
	      if(ncmpl==1) then
	        call boson_unpack_tau3r(tmm,tmi,tii,nplw_gw,
     &                              p_w_q_tau_mm(1,1,ind_k),
     &                              p_w_q_tau_mi(1,1,1,ind_k),
     &                              p_w_q_tau_ii(1,1,ind_k),nplw_gw,
     &                              nplwgw(k),k)
          else
            call boson_unpack_tau3(tmm,tmi,tii,nplw_gw,
     &                             p_w_q_tau_mm(1,1,ind_k),
     &                             p_w_q_tau_mi(1,1,1,ind_k),
     &                             p_w_q_tau_ii(1,1,ind_k),nplw_gw,
     &                             nplwgw(k))
          endif
	      call check_w_sym(k,tmm,tmi,tii,nplw_gw)
        enddo
      endif
      call ws_k_nu_solid
	  if(itest/=0.and..not.lsclp) then
	    if(maswrk) write(iun,*)' Symmetry check for W:'
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n_pbt=n_pbmt+nplwgw(k)
	      if(ncmpl==1) then
	        call boson_unpack_tau3r(tmm,tmi,tii,nplw_gw,
     &                              p_w_q_tau_mm(1,1,ind_k),
     &                              p_w_q_tau_mi(1,1,1,ind_k),
     &                              p_w_q_tau_ii(1,1,ind_k),nplw_gw,
     &                              nplwgw(k),k)
          else
            call boson_unpack_tau3(tmm,tmi,tii,nplw_gw,
     &                             p_w_q_tau_mm(1,1,ind_k),
     &                             p_w_q_tau_mi(1,1,1,ind_k),
     &                             p_w_q_tau_ii(1,1,ind_k),nplw_gw,
     &                             nplwgw(k))
          endif
	      call check_w_sym(k,tmm,tmi,tii,nplw_gw)
        enddo
      endif
      if(itest/=0) deallocate(tmm,tmi,tii)
      call timel('**** SCREENED_W finished ***********')
      end
