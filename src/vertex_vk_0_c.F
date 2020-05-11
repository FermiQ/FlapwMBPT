      subroutine vertex_vk_0_c(wpb,k_rs,vx,fif,n,np,n0)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,np,n0
      real*8, intent(in) :: wpb(n0,n0)
      complex*16, intent(in) :: k_rs(n,n),fif(n,n,np)
      complex*16, intent(inout) :: vx(n,n)
      integer :: i,j,ii,jj
      complex*16, allocatable :: aa(:,:,:),a0(:,:,:,:)
      allocate(aa(n,n,np))
      call dgemm('n','n',2*n*n,np,np,1.d0,fif,2*n*n,wpb,n0,0.d0,aa,
     &           2*n*n)
      allocate(a0(n,n,n,n))
      call zgemm('n','c',n*n,n*n,np,(1.d0,0.d0),aa,n*n,fif,n*n,
     &           (0.d0,0.d0),a0,n*n)
      vx=(0.d0,0.d0)
      do j=1,n
        do i=1,n
          do jj=1,n
            do ii=1,n
              vx(i,j)=vx(i,j)+k_rs(ii,jj)*a0(i,ii,j,jj)
            enddo
          enddo
        enddo
      enddo
      deallocate(aa,a0)
      end
