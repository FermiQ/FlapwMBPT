      subroutine v_loc_from_red(vloc,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: khbr
      real*8, intent(out) :: vloc(n_pbmtm_red,n_pbmtm_red,natom)
      integer :: k,ig,ind_q,iq,n,kst,nn,iatom,i0,isort,j,i,ka
      real*8 :: qq
      complex*16, allocatable :: tmp(:,:)
      n=n_pbtot_red
      nn=n*n
      qq=8.d0*pi*q2aver_c
      if(key_hybrid==1) qq=qq*alpha_hbr
c ----------------------------------------------------------------
      allocate(tmp(n,n))
      vloc=0.d0
      do ind_q=1,ndim_k_red(me_k+1)
        k=n_mpi_k_red(me_k+1)+ind_q
        ka=k_a_from_c(k)
        do kst=1,k_star_c(k)
          iq=k_list_c(kst,k)
          ig=k_group_c(iq)
          if(khbr==0) tmp=v_red_q(:,:,ka)
          if(khbr==1) tmp=v_red_b(:,:,k)
          call sym_w_red(iq,k,tmp)
          do iatom=1,natom
            i0=iopb_red(iatom)-1
            isort=is(iatom)
            do j=1,n_pbmt0_red(isort)
              do i=1,n_pbmt0_red(isort)
                vloc(i,j,iatom)=vloc(i,j,iatom)+tmp(i0+i,i0+j)
                if(k==1.and.khbr==0) then
c --------- Add singular q=0 term ----------------------------------
                  vloc(i,j,iatom)=vloc(i,j,iatom)+pw_pb_red(i0+i,1)*qq
     &                                        *conjg(pw_pb_red(i0+j,1))
                endif
              enddo
            enddo
          enddo
        enddo   !! over kst
      enddo   !! over ind_q
      nn=n_pbmtm_red**2
      deallocate(tmp)
      if(nproc_k/=1) call dgop(vloc,nn*natom,'  +',comm_k)
      vloc=vloc/nqdiv_c
      end
