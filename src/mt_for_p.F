      subroutine mt_for_p(tt,npw,nsta,istar,ind_k)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nsta,istar,npw,ind_k
      complex*16, intent(out) :: tt(nsta,npw)
      integer :: i,i0,j,j0,k
      real*8 :: v(3),phase,v1(3),pi2
      real*8, allocatable :: sr(:,:)
      complex*16, allocatable :: p(:,:),s(:,:)
      k=n3_mpi_k(me_k+1)+ind_k
      pi2=pi+pi
      allocate(p(nsta,npw))
      do i0=1,nsta
        i=ind_r_star(i0,istar)
        v=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)+r_pnt(3,i)*r0b0(:,3)
        do j=1,npw
          j0=indpw_gw(j,k)
          v1=pnt(:,k)+gbs(:,j0)
          phase=pi2*dot_product(v1,v)
          p(i0,j)=dcmplx(cos(phase),sin(phase))
        enddo
      enddo
      p=p*sqrt(amega)/nr_full
      if(complex_ro) then
        allocate(s(npw,npw))
        call s_overlap(0,s,npw,npw,indpw_gw(1,k))
        s=s/amega
        call zgemm('n','n',nsta,npw,npw,(1.d0,0.d0),p,nsta,s,npw,
     &           (0.d0,0.d0),tt,nsta)
        deallocate(s)
      else
        allocate(sr(npw,npw))
        call s_overlap_r(0,sr,npw,npw,indpw_gw(1,k))
        call dgemm('n','n',2*nsta,npw,npw,1.d0,p,2*nsta,sr,npw,0.d0,tt,
     &             2*nsta)
        deallocate(sr)
      endif
      deallocate(p)
      end
