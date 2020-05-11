      subroutine w4_loc_stat(w4)
      use atom_mod
      use heg_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(out) :: w4(maxel_red**2,maxel_red**2,natom,nspin,
     &                          nspin)
      integer :: k,ig,ind_q,iq,n,kst,nn,iatom,i0,isort,j,i,ind_nu,i_nu,
     &           mx2,ispin,ij,i1,n1,n11,jspin,ka
      real*8 :: ww
      real*8, allocatable :: vloc(:,:,:),fifim(:,:,:,:),fif(:,:,:,:),
     &                       tmp33(:,:)
	  complex*16, allocatable :: tmp(:,:)
      mx2=maxel_red**2
      n=n_pbtot_red
      allocate(tmp(n,n))
      allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
      vloc=0.d0
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
                  vloc(i,j,iatom)=vloc(i,j,iatom)+tmp(i0+i,i0+j)
                  if(k==1) then
c --------- Add singular q=0 term ----------------------------------
                    vloc(i,j,iatom)=vloc(i,j,iatom)
     &                                  +pw_pb_tild_red(i0+i,1)*ww
     &                                  *conjg(pw_pb_tild_red(i0+j,1))
                  endif
                enddo
              enddo
            enddo
          enddo   !! over kst
        enddo   !! over ind_q
      enddo
      deallocate(tmp)
      nn=n_pbmtm_red**2
	  if(nproc_t/=1) call dgop(vloc,nn*natom,'  +',comm_t)
	  if(nproc_k/=1) call dgop(vloc,nn*natom,'  +',comm_k)
      vloc=vloc/nqdiv_c
c ---------------------------------------------------------------
      allocate(fif(mx2,n_pbmtm_red,nsort,nspin))
      allocate(fifim(maxel_red,maxel_red,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do isort=1,nsort
          ij=0
          do i1=1,lfunm_red(isort)
            do i=1,lfunm_red(isort)
              ij=ij+1
              fif(ij,:,isort,ispin)=fifim(i,i1,:,isort)
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
      allocate(tmp33(mx2,n_pbmtm_red))
	  do iatom=1,natom
	    isort=is(iatom)
	    n1=lfunm_red(isort)
	    n11=n_pbmt0_red(isort)
        do jspin=1,nspin
	      do ispin=1,nspin
            call dgemm('n','n',n1*n1,n11,n11,1.d0,fif(1,1,isort,ispin),
     &                 mx2,vloc(1,1,iatom),n_pbmtm_red,0.d0,tmp33,mx2)
	        call dgemm('n','t',n1*n1,n1*n1,n11,1.d0,tmp33,mx2,
     &	               fif(1,1,isort,jspin),mx2,0.d0,
     &                 w4(1,1,iatom,ispin,jspin),mx2)
          enddo   !! over ispin
        enddo   !! over jspin
      enddo  !! over iatom
      deallocate(vloc,tmp33,fif)
      end
