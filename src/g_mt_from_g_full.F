      subroutine g_mt_from_g_full(correlation,gloc)
	  use atom_mod
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      logical, intent(in) :: correlation
      complex*16 :: gloc(maxel,maxel,natom,nspin)
	  integer :: k,iatom,isort,ie,ind_k,n,ispin,i,ii
	  real*8 :: de,green0,ee
      complex*16 :: wei
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:)
	  gloc=(0.d0,0.d0)
      allocate(tmp(maxel,nbndf))
      allocate(tmp1(nbndf,nbndf))
      allocate(tmp2(nbndf,nbndf))
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
	      wei=dcmplx(wgt(k),0.d0)
	      tmp1=(0.d0,0.d0)
	      if(.not.correlation) then
            if(k_integral=='FD') then
              do i=1,n
                de=e_bnd(i,k,ispin)-chem_pot
                tmp1(i,i)=-green0(de,betta_t)
              enddo
            else if(k_integral=='O1') then
              do ii=1,nn_k_int
                do i=1,n
                  ee=e_bnd(i,k,ispin)
     &              +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &              +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &              +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
                  de=ee-chem_pot
                  tmp1(i,i)=tmp1(i,i)-green0(de,betta_t)
                enddo
              enddo
              tmp1=tmp1/nn_k_int
            endif
	      else
            if(me_t==0) then
              call ferm_unpack_tau(tmp1,g_full(1,1,1,1,ind_k,ispin),n,
     &                             nbndf,nbndf,2)
              tmp1=-tmp1
            endif
            if(nproc_t/=1) call dgop(tmp1,2*nbndf**2,'  +',comm_t)
          endif
	      do iatom=1,natom
	        isort=is(iatom)
	        ie=io_lem(iatom)
	        call zgemm('n','n',lfunm(isort),n,n,(1.d0,0.d0),
     &                 z_bnd(ie,1,ind_k,ispin),nfun,tmp1,nbndf,
     &                 (0.d0,0.d0),tmp,maxel)
	        call zgemm('n','c',lfunm(isort),lfunm(isort),n,wei,tmp,
     &                 maxel,z_bnd(ie,1,ind_k,ispin),nfun,(1.d0,0.d0),
     &                 gloc(1,1,iatom,ispin),maxel)
	      enddo   !! over iatom
	    enddo   !! over ind_k
	  enddo   !! over ispin
	  deallocate(tmp,tmp1,tmp2)
      if(nproc_k/=1) call DGOP(gloc,2*maxel**2*natom*nspin,'  +',comm_k)
c ------- Symmetrization -----------------------------------------------
      do ispin=1,nspin
	    call symg_loc(gloc(1,1,1,ispin))
	    gloc(:,:,:,ispin)=conjg(gloc(:,:,:,ispin))
	  enddo
      end
