      subroutine tsig_loc_c(gtau,k_rs,t_rs,fif,n,np)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,np
      complex*16, intent(in) :: fif(n,n,np),gtau(n,n),k_rs(n,np)
      complex*16, intent(out) :: t_rs(np,n)
      integer :: i,j,nn,ii
      complex*16 :: tt
      complex*16, allocatable :: aa(:,:,:),bb(:,:,:),dd(:,:)
      allocate(aa(n,n,n))
      aa=(0.d0,0.d0)
      do j=1,n
        do i=1,n
          do ii=1,np
            tt=fif(i,j,ii)
            aa(:,i,j)=aa(:,i,j)+k_rs(:,ii)*tt
          enddo
        enddo
      enddo
      allocate(bb(n,n,n))
      nn=n*n
      call zgemm('n','c',nn,n,n,(1.d0,0.d0),aa,nn,gtau,n,(0.d0,0.d0),bb,
     &           nn)
      deallocate(aa)
      allocate(dd(n,np))
      dd=(0.d0,0.d0)
      do j=1,n
        do i=1,n
          do ii=1,np
            tt=fif(i,j,ii)
            dd(:,ii)=dd(:,ii)+bb(j,:,i)*tt
          enddo
        enddo
      enddo
      do j=1,n
        do i=1,np
          t_rs(i,j)=-dd(j,i)
        enddo
      enddo
      deallocate(bb,dd)
      end
