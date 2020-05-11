      subroutine chi0_lda_uniform(isp1,isp2,chi0,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
#ifdef MPI
	include 'mpif.h'
#endif
      integer, intent(in) :: isp1,isp2,n
      complex*16, intent(out) :: chi0(n,n)
      integer :: k,ind_k,i,j,jj,ispin,ig,n1,n2,npr,me1,n12,j0
      real*8 :: gtild(3,2),de,green0,aa
      integer, allocatable :: n2d(:),n20(:)
      real*8, allocatable :: gb(:,:),fifim(:,:,:,:)
      complex*16, allocatable :: ps3(:,:,:),tmp(:,:),tmp1(:,:),t(:,:,:)
      npr=nproc_t*nproc_b
      me1=me_t_b+1
      allocate(n2d(npr))
      allocate(n20(npr))
      allocate(tmp(n,n))
      tmp=(0.d0,0.d0)
      gtild=0.d0            
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
      call fi0_full(fifim,isp1,isp2)
      allocate(gb(nbndf,nspin))
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
        do ispin=1,nspin
          do i=1,n_bnd(k,ispin)
            de=e_bnd(i,k,ispin)-chem_pot
            gb(i,ispin)=green0(de,betta_t)
          enddo
        enddo
        n1=n_bnd(k,isp1)
        n2=n_bnd(k,isp2)
        call size_shift_par(n2,npr,n2d,n20)
        allocate(ps3(n1,n2d(me1),n))
        call ppm_factors0(k,k,1,fifim,z_bnd(1,1,ind_k,isp1),
     &	                  z_bnd(1,1,ind_k,isp2),
     &                    ev_bnd(1,1,ind_k,isp1),
     &                    ev_bnd(1,1,ind_k,isp2),ps3,gtild,n,n1,
     &                    n2d(me1),n20(me1),nbndf)
        allocate(t(n1,n2d(me1),n))
        do j=1,n2d(me1)
          j0=n20(me1)+j  
          do i=1,n1
            de=e_bnd(j0,k,isp2)-e_bnd(i,k,isp1)
            if(abs(de)<1.d-8) then
              aa=betta_t*(1.d0+gb(i,isp1))*gb(i,isp1)
            else
              aa=(gb(i,isp1)-gb(j0,isp2))/de
            endif
            do jj=1,n
              t(i,j,jj)=aa*ps3(i,j,jj)
            enddo
          enddo
        enddo
        t=wgt(k)*t
        n12=n1*n2d(me1)
        call zgemm('c','n',n,n,n12,(1.d0,0.d0),ps3,n12,t,n12,
     &             (1.d0,0.d0),tmp,n)
        deallocate(ps3,t)
      enddo   !! over ind_k
      deallocate(fifim,gb,n2d,n20)
      if(nproc_k/=1) call dgop(tmp,2*n*n,'  +',comm_k)
      if(npr/=1) call dgop(tmp,2*n*n,'  +',comm_t_b)
c --------- Symmetrization --------------------------------------
      chi0=(0.d0,0.d0)
      allocate(tmp1(n,n))
      allocate(n2d(nproc))
      allocate(n20(nproc))
      call size_shift_par(ngroup,nproc,n2d,n20)
      me1=me+1
      do i=1,n2d(me1)
        ig=n20(me1)+i
        tmp1=tmp
        call sym_chi_gamma(ig,tmp1,n)
        chi0=chi0+tmp1
      enddo
      chi0=chi0/ngroup
      deallocate(tmp,tmp1,n2d,n20)
#ifdef MPI
      if(nproc/=1) call dgop(chi0,2*n*n,'  +',MPI_COMM_WORLD)
#endif	  
      call timel('**** CHI0_LDA_UNIFORM finished *****')
      end
