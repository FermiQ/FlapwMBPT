      subroutine sig_psi5_g2_w1_stat(tmatr,sig,g_rr_tau,ch)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      real*8, intent(in) :: g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,
     &                               ndim3_tau,nspin)
      complex*16, intent(in) :: tmatr(n_lem2_red,n_lem2_red,ndim3_nu,
     &                                ndimc_kk(me_k+1),nspin,nspin)
      complex*16, intent(out) :: sig(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                               ndim_k_red(me_k+1),nspin)
      integer :: iatom,isort,ind2,n,nn,ispin,ir,ind_r,ind_tau,it,
     &           jatom,jsort,n1,nn1,jnd,jnd2,ind,l,i,ind_nu,ind_q,iq,
     &           iq1,jspin,k,ind_k,ka,kl,ij,j
      real*8 :: pi2,phase
      complex*16 :: cc
      real*8, allocatable :: t_rr(:,:,:,:,:,:,:),sigr(:,:,:,:,:,:)
      complex*16, allocatable :: kr(:,:,:),tn(:,:,:,:),sigk(:,:,:),
     &                           tmp(:,:)
      pi2=pi+pi
c ---   T-matrix in (R;tau) veriables ------------------------------
      allocate(t_rr(n_lem2_red,n_lem2_red,2,ndim3_tau,
     &              ndimc_kk(me_k+1),nspin,nspin))
      allocate(tn(n_lem2_red,n_lem2_red,ndim3_nu,ndimc_kk(me_k+1)))
	  allocate(kr(nqdiv_c,n_lem2_red,n_lem2_red))
	  do jspin=1,nspin
        do ispin=1,nspin
          do ind_nu=1,ndim3_nu
            kr=(0.d0,0.d0)
	        do ind_q=1,ndimc_kk(me_k+1)
	          iq=nc_mpi_kk(me_k+1)+ind_q
              call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,iq1)
              kr(iq1,:,:)=tmatr(:,:,ind_nu,ind_q,ispin,jspin)
            enddo
            if(nproc_k/=1) call dgop(kr,2*n_lem2_red**2*nqdiv_c,'  +',
     &                               comm_k)
            call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_lem2_red**2,kr,1)
	        do ind_r=1,ndimc_kk(me_k+1)
	          ir=nc_mpi_kk(me_k+1)+ind_r
	          tn(:,:,ind_nu,ind_r)=kr(ir,:,:)/nqdiv_c
	        enddo
	      enddo   !! over ind_nu
	      do ind_r=1,ndimc_kk(me_k+1)
	        call nu_to_tau_c_spl(tn(1,1,1,ind_r),
     &	                         t_rr(1,1,1,1,ind_r,ispin,jspin),
     &	                         n_lem2_red**2)
          enddo
        enddo
      enddo      
      if(nspin==1.and.irel/=2) t_rr=2.d0*t_rr
      deallocate(kr,tn)
c ----------------------------------------------------------------
      allocate(sigr(nfun_red,nfun_red,2,ndim3_tau,ndimc_kk(me_k+1),
     &              nspin))
      sigr=0.d0
	  do jatom=1,natom
	    jsort=is(jatom)
	    n1=lfunm_red(jsort)
	    nn1=n1*n1
        jnd=io_lem_red(jatom)-1
        jnd2=io_lem2_red(jatom)-1
	    do iatom=1,natom
	      isort=is(iatom)
	      n=lfunm_red(isort)
	      nn=n*n
          ind=io_lem_red(iatom)-1
          ind2=io_lem2_red(iatom)-1
	      do jspin=1,nspin
            do ispin=1,nspin
              do ind_tau=1,ndim3_tau
                do it=1,2
                  do ind_r=1,ndimc_kk(me_k+1)
	                ir=nc_mpi_kk(me_k+1)+ind_r
                    kl=0
                    do l=1,n1
                      do k=1,n1
                        kl=kl+1
                        ij=0
                        do j=1,n
                          do i=1,n
                            ij=ij+1
                            if(ch=='ph') then
                              sigr(ind+i,jnd+k,it,ind_tau,ind_r,ispin)=
     &	                        sigr(ind+i,jnd+k,it,ind_tau,ind_r,ispin)
     &                         +g_rr_tau(ind+j,jnd+l,ir,it,ind_tau,
     &                                   jspin)
     &                         *t_rr(ind2+ij,jnd2+kl,it,ind_tau,ind_r,
     &                               ispin,jspin)
                            else
                              sigr(ind+i,jnd+k,it,ind_tau,ind_r,ispin)=
     &	                        sigr(ind+i,jnd+k,it,ind_tau,ind_r,ispin)
     &                         -g_rr_tau(ind+j,jnd+l,ir,3-it,ind_tau,
     &                                   jspin)
     &                         *t_rr(ind2+ij,jnd2+kl,it,ind_tau,ind_r,
     &                               ispin,jspin)
                            endif
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo  !! over ind_r
                enddo  !! over it
              enddo  !! over ind_tau
            enddo   !! over ispin
          enddo  !! over jspin
        enddo  !! over iatom_c
      enddo  !! over jatom_c
      deallocate(t_rr)
      allocate(sigk(nfun_red,nfun_red,npnt_c))
      allocate(tmp(nfun_red,nbndf_bnd))
      do ispin=1,nspin
        do ind_tau=1,ndim3_tau
          do it=1,2
            sigk=(0.d0,0.d0)
            do ind_r=1,ndimc_kk(me_k+1)
	        ir=nc_mpi_kk(me_k+1)+ind_r
              do k=1,npnt_c
                phase=-pi2*dot_product(pnt_c(:,k),rr_coarse(:,ir))
                cc=dcmplx(cos(phase),sin(phase))
                sigk(:,:,k)=sigk(:,:,k)
     &                     +cc*sigr(:,:,it,ind_tau,ind_r,ispin)
              enddo
            enddo
            if(nproc_k/=1) call dgop(sigk,2*nfun_red**2*npnt_c,'  +',
     &                               comm_k)
	        do ind_k=1,ndim_k_red(me_k+1)
	          k=n_mpi_k_red(me_k+1)+ind_k
              ka=k_a_from_c(k)
              n=n_low_bnd(ka,ispin)
	          call zgemm('n','n',nfun_red,n,nfun_red,(1.d0,0.d0),
     &                   sigk(1,1,k),nfun_red,z_bnd_red(1,1,ka,ispin),
     &                   nfun_red,(0.d0,0.d0),tmp,nfun_red)
              call zgemm('c','n',n,n,nfun_red,(1.d0,0.d0),
     &                   z_bnd_red(1,1,ka,ispin),nfun_red,tmp,
     &	                 nfun_red,(0.d0,0.d0),
     &                   sig(1,1,it,ind_tau,ind_k,ispin),nbndf_bnd)
            enddo   !! over ind_k
          enddo  !! over it
        enddo   !! over ind_tau
      enddo   !! over ispin
      deallocate(sigr,sigk,tmp)
      call timel('* SIG_PSI5_G2_W1 finished **********')
      end
