      subroutine a_z_symmetrize(ind_k,abnd,zbnd,nbn)
      use atom_mod
	use manager_mod 
	use parallel_mod    
	use solid_mod
      implicit none
      integer, intent(in) :: ind_k,nbn
      complex*16, intent(inout) :: abnd(nbasmpw,nbndf),zbnd(nfun,nbndf)
      integer :: k,ig,ig0
      complex*16, allocatable :: aold(:,:),zold(:,:),an(:,:),zn(:,:)
      k=n3_mpi_k(me_k+1)+ind_k
      allocate(aold(nbasmpw,nbn))
      allocate(zold(nfun,nbn))
      allocate(an(nbasmpw,nbn))
      allocate(zn(nfun,nbn))
      aold=abnd(:,1:nbn)
      zold=zbnd(:,1:nbn)
      abnd=(0.d0,0.d0)
      zbnd=(0.d0,0.d0)
	do ig0=1,num_sym_k(k)
        ig=list_sym_k(ig0,k)
        call sym_z(zn,zold,ig,pnt(1,k),nfun,nbn,maxb,io_lem,lmb,indbasa,
     &             limlb,0)
        call sym_a(an,aold,k,ig,pnt(1,k),nbasmpw,nbn,nbask(k),
     &             indgb(1,k),iplf_bk(1,k))
        abnd(:,1:nbn)=abnd(:,1:nbn)+an
        zbnd(:,1:nbn)=zbnd(:,1:nbn)+zn
      enddo
      abnd=abnd/num_sym_k(k)
      zbnd=zbnd/num_sym_k(k)
      deallocate(aold,zold,an,zn)
      end