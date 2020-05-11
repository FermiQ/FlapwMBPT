      subroutine l_q_nu_g2_w1(ispin,jspin,l_pp,g_rr_tau,ch)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: ispin,jspin
      real*8, intent(in) :: g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,
     &                               ndim3_tau,nspin)
      complex*16, intent(out) :: l_pp(n_lem2_red,n_lem2_red,
     &                                ndim3_nu,ndimc_kk(me_k+1))
	  integer :: ind_tau,it,jatom,jsort,n1,jnd,iatom,isort,ind,n,ir,l,
     &           jnd2,ind2,ind_r,ind_nu,ind_q,iq,i,j,k,ik,jl
      real*8, allocatable :: l_rr_t(:,:,:,:,:)
	  complex*16, allocatable :: kr(:,:,:)
	  allocate(l_rr_t(n_lem2_red,n_lem2_red,2,ndim3_tau,
     &                ndimc_kk(me_k+1)))
	  do jatom=1,natom
	    jsort=is(jatom)
	    n1=lfunm_red(jsort)
        jnd=io_lem_red(jatom)-1
        jnd2=io_lem2_red(jatom)-1
	    do iatom=1,natom
	      isort=is(iatom)
	      n=lfunm_red(isort)
          ind=io_lem_red(iatom)-1
          ind2=io_lem2_red(iatom)-1
          do ind_tau=1,ndim3_tau
            do it=1,2
              do ind_r=1,ndimc_kk(me_k+1)
	            ir=nc_mpi_kk(me_k+1)+ind_r
                jl=0
	            do l=1,n1
	              do j=1,n1
                    jl=jl+1
                    ik=0
	                do k=1,n
	                  do i=1,n
                        ik=ik+1
	                    if(ch=='pp') then
	                      l_rr_t(ind2+ik,jnd2+jl,it,ind_tau,ind_r)=
     & 	                    g_rr_tau(ind+i,jnd+j,ir,it,ind_tau,ispin)
     &	                   *g_rr_tau(ind+k,jnd+l,ir,it,ind_tau,jspin)
	                    else if(ch=='ph') then
	                      l_rr_t(ind2+ik,jnd2+jl,it,ind_tau,ind_r)=
     & 	                   -g_rr_tau(ind+i,jnd+j,ir,it,ind_tau,ispin)
     &	                   *g_rr_tau(ind+k,jnd+l,ir,3-it,ind_tau,
     &                                 jspin)
                        endif
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c ------------------------------------------------------------------
      do ind_r=1,ndimc_kk(me_k+1)
        call tau_to_nu_re(l_pp(1,1,1,ind_r),l_rr_t(1,1,1,1,ind_r),
     &                    n_lem2_red**2)
      enddo
      deallocate(l_rr_t)
c ------------------------------------------------------------------
	  allocate(kr(nqdiv_c,n_lem2_red,n_lem2_red))
      do ind_nu=1,ndim3_nu
        kr=(0.d0,0.d0)
        do ind_r=1,ndimc_kk(me_k+1)
	      ir=nc_mpi_kk(me_k+1)+ind_r
	      kr(ir,:,:)=l_pp(:,:,ind_nu,ind_r)
	    enddo
        if(nproc_k/=1) call dgop(kr,2*n_lem2_red**2*nqdiv_c,'  +',
     &                           comm_k)
	    call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_lem2_red**2,kr,-1)
        do ind_q=1,ndimc_kk(me_k+1)
	      iq=nc_mpi_kk(me_k+1)+ind_q
          call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,ir)
          l_pp(:,:,ind_nu,ind_q)=kr(ir,:,:)
        enddo
      enddo
	  deallocate(kr)
	  end
