      subroutine w_loc_from_red(vloc,wloc)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(out) :: wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu,natom),
     &                       vloc(n_pbmtm_red,n_pbmtm_red,natom)
      integer :: k,ig,ind_q,iq,n,kst,ind_nu,i_nu,nn,iatom,i0,isort,j,i,
     &           ka
	  real*8 :: qq,ww
	  complex*16, allocatable :: tmp(:,:)
      n=n_pbtot_red
      nn=n*n
      qq=8.d0*pi*q2aver_c
c ----------------------------------------------------------------
      allocate(tmp(n,n))
      wloc=0.d0
      vloc=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        k=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(k)
	    do kst=1,k_star_c(k)
	      iq=k_list_c(kst,k)
	      ig=k_group_c(iq)
          tmp=v_red_q(:,:,ka)
	      call sym_w_red(iq,k,tmp)
          do iatom=1,natom
            i0=iopb_red(iatom)-1
            isort=is(iatom)
            do j=1,n_pbmt0_red(isort)
              do i=1,n_pbmt0_red(isort)
                vloc(i,j,iatom)=vloc(i,j,iatom)+tmp(i0+i,i0+j)
                if(k==1) then
c --------- Add singular q=0 term ----------------------------------
                  vloc(i,j,iatom)=vloc(i,j,iatom)
     &                           +pw_pb_tild_red(i0+i,1)*qq
     &                           *conjg(pw_pb_tild_red(i0+j,1))
                endif
              enddo
            enddo
          enddo
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            tmp=w_red_q(:,:,ind_nu,ka)
	        call sym_w_red(iq,k,tmp)
	        if(k==1) ww=ws_head_nu(i_nu)*q2aver_c
            do iatom=1,natom
              i0=iopb_red(iatom)-1
              isort=is(iatom)
              do j=1,n_pbmt0_red(isort)
                do i=1,n_pbmt0_red(isort)
                  wloc(i,j,i_nu,iatom)=wloc(i,j,i_nu,iatom)
     &                                +tmp(i0+i,i0+j)
                  if(k==1) then
c --------- Add singular q=0 term ----------------------------------
                    wloc(i,j,i_nu,iatom)=wloc(i,j,i_nu,iatom)
     &                                  +pw_pb_tild_red(i0+i,1)*ww
     &                                  *conjg(pw_pb_tild_red(i0+j,1))
                  endif
                enddo
              enddo
            enddo
          enddo  !! over ind_nu
        enddo   !! over kst
      enddo   !! over ind_q
      nn=n_pbmtm_red**2
	  if(nproc_t/=1) call dgop(wloc,nn*(n_nu+1)*natom,'  +',
     & 	                         comm_t)
      deallocate(tmp)
	  if(nproc_k/=1) then
	    call dgop(wloc,nn*(n_nu+1)*natom,'  +',comm_k)
	    call dgop(vloc,nn*natom,'  +',comm_k)
      endif
      wloc=wloc/nqdiv_c
      vloc=vloc/nqdiv_c
      end
