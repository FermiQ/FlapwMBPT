      subroutine bands_qp(x)
c     x - current chemical potential
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: x
      integer :: ispin,k,n,j,ind_k,i,il
      real*8, allocatable :: sigc(:,:,:,:)
	  complex*16, allocatable :: a(:,:),b(:,:),cc(:,:)
c ------------------------------------------------------------------
      e_qp=0.d0
      z_qp=(0.d0,0.d0)
	  do ispin=1,nspin
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n=n_bnd(k,ispin)
          allocate(sigc(n,n,2,ndim3_omega))
          call tau_to_omega_baa(sigc,n,sig_c_tau(1,1,1,1,ind_k,ispin),
     &                          nbndf,n)
          allocate(a(n,n))
          allocate(b(n,n))
          call sig_omega0_solid(n,sigc,ind_k,ispin,a,b)
          deallocate(sigc)
c  Matrix b is inverse Z factor now - we diagonalize it to get Z^(1/2)
c ---------- Eigen pairs of Z^(-1) ----------------------------------
          allocate(cc(n,n))
	      do j=1,n
	        do i=1,n
	          cc(i,j)=b(i,j)
	        enddo
	      enddo
          if(me_t_b==0) call eig_val_solver(n,n,cc,e_qp(1,k,ispin))
          if(nproc_t*nproc_b/=1) then
	        call brdcst(msgdbl,e_qp(1,k,ispin),8*n,0,comm_t_b)
	        call brdcst(msgdbl,cc,16*n*n,0,comm_t_b)
          endif
          if(e_qp(1,k,ispin)<=0.d0) then
            print *,'me : ',me,' problem in BANDS_QP'
            print *,' k diag(1)',k,e_qp(1,k,ispin)
          endif
          if(nbndf>n) e_qp(n+1:nbndf,k,ispin)=e_qp(n,k,ispin)+0.2d0
	      b=(0.d0,0.d0)
	      do j=1,n
	        do i=1,n
	          do il=1,n
	            b(i,j)=b(i,j)+cc(i,il)/sqrt(e_qp(il,k,ispin))*
     &                          conjg(cc(j,il))
	          enddo
	        enddo
	      enddo
c ----- Now B = Z^(1/2) ------------------------------------------     
	      do j=1,n
	        a(j,j)=a(j,j)+e_bnd(j,k,ispin)-x
	      enddo
          call zgemm('n','n',n,n,n,(1.d0,0.d0),a,n,b,n,(0.d0,0.d0),cc,n)
          call zgemm('n','n',n,n,n,(1.d0,0.d0),b,n,cc,n,(0.d0,0.d0),a,n)
	      deallocate(b)
	      do j=1,n
	        do i=1,n
	          cc(i,j)=a(i,j)
	        enddo
	      enddo
	      deallocate(a)
          if(me_t_b==0) call eig_val_solver(n,n,cc,e_qp(1,k,ispin))
          if(nproc_t*nproc_b/=1) then
	        call brdcst(msgdbl,e_qp(1,k,ispin),8*n,0,comm_t_b)
	        call brdcst(msgdbl,cc,16*n*n,0,comm_t_b)
          endif
	      do j=1,n
	        e_qp(j,k,ispin)=e_qp(j,k,ispin)+x
	      enddo
c --- We put renormalization Z = 1 here ------------------------	    
	      q_qp(1:n,1:n,ind_k,ispin)=cc
	      deallocate(cc)
c ---------------------------------------------------------------	    
          call zgemm('n','n',nfun,n,n,(1.d0,0.d0),
     &               z_bnd(1,1,ind_k,ispin),nfun,q_qp(1,1,ind_k,ispin),
     &               nbndf,(0.d0,0.d0),z_qp(1,1,ind_k,ispin),nfun)
          do i=1,n
            qp_spectr_weight(i,ind_k,ispin)=
     &          sum(abs(q_qp(1:n,i,ind_k,ispin))**2)
          enddo
	    enddo   ! over ind_k
	  enddo  !! over ispin
	  if(nproc_k/=1) call dgop(e_qp,nbndf*npnt*nspin,'  +',comm_k)
      end
