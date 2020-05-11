      subroutine w2_loc_stat(w2)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(out) :: w2(n_pbmtm_red,n_pbmtm_red,natom)
      integer :: k,ig,ind_q,iq,n,kst,nn,iatom,i0,isort,j,i,ind_nu,i_nu,
     &           ka
      real*8 :: ww
      complex*16, allocatable :: tmp(:,:)
      n=n_pbtot_red
      allocate(tmp(n,n))
      w2=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(i_nu/=nu_w_stat) cycle
        do ind_q=1,ndim_k_red(me_k+1)
          k=n_mpi_k_red(me_k+1)+ind_q
          if(k==1) ww=(8.d0*pi+ws_head_nu(i_nu))*q2aver_c
          ka=k_a_from_c(k)
          do kst=1,k_star_c(k)
            iq=k_list_c(kst,k)
            ig=k_group_c(iq)
            tmp=v_red_q(:,:,ka)+w_red_q(:,:,ind_nu,ka)
            call sym_w_red(iq,k,tmp)
            do iatom=1,natom
              i0=iopb_red(iatom)-1
              isort=is(iatom)
              do j=1,n_pbmt0_red(isort)
                do i=1,n_pbmt0_red(isort)
                  w2(i,j,iatom)=w2(i,j,iatom)+tmp(i0+i,i0+j)
                  if(k==1) then
c --------- Add singular q=0 term ----------------------------------
                    w2(i,j,iatom)=w2(i,j,iatom)+pw_pb_tild_red(i0+i,1)
     &                           *ww*conjg(pw_pb_tild_red(i0+j,1))
                  endif
                enddo
              enddo
            enddo
          enddo   !! over kst
        enddo   !! over ind_q
      enddo
      deallocate(tmp)
      nn=n_pbmtm_red**2
      if(nproc_t/=1) call dgop(w2,nn*natom,'  +',comm_t)
      if(nproc_k/=1) call dgop(w2,nn*natom,'  +',comm_k)
      w2=w2/nqdiv_c
      end
