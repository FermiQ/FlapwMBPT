      subroutine vertex_loc_gt_c(i_tau,ll,t_rs,gtau,vx,fif,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll,n,np
      complex*16, intent(in) :: gtau(n,n,0:n_tau),fif(n,n,np),
     &                          t_rs(np,np,ll)
      complex*16, intent(inout) :: vx(n,n,ll)
      integer :: i,j,il,j1,i1,il3
      complex*16, allocatable :: cc(:,:),aa(:,:,:),bb(:,:,:,:)
      allocate(aa(np,n,n))
      allocate(bb(n,n,n,n))
      allocate(cc(n,n))
      do il=1,ll
        il3=3-il
        if(ll==1) il3=1
        call zgemm('n','c',np,n*n,np,(1.d0,0.d0),t_rs(1,1,il3),np,fif,
     &             n*n,(0.d0,0.d0),aa,np)
        call zgemm('n','n',n*n,n*n,np,(1.d0,0.d0),fif,n*n,aa,np,
     &             (0.d0,0.d0),bb,n*n)
        cc=(0.d0,0.d0)
        do j1=1,n
          do i1=1,n
            do j=1,n
              do i=1,n
                cc(i,j)=cc(i,j)+gtau(i1,j1,i_tau)*bb(j1,j,i1,i)
              enddo
            enddo
          enddo
        enddo
        vx(:,:,il)=vx(:,:,il)+cc
      enddo   !! over il
      deallocate(aa,bb,cc)
      end
