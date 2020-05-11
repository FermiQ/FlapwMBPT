      subroutine vertex_loc_gt_r(i_tau,ll,t_rs,gtau,vx,fif,n,np)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll,n,np
      real*8, intent(in) :: gtau(n,n,0:n_tau),fif(n,n,np)
      complex*16, intent(in) :: t_rs(np,np,ll)
      complex*16, intent(inout) :: vx(n,n,ll)
      integer :: i,j,il,j1,i1,il3,it
      real*8, allocatable :: cc(:,:),aa(:,:,:),bb(:,:,:,:),tt(:,:)
      allocate(aa(np,n,n))
      allocate(bb(n,n,n,n))
      allocate(cc(n,n))
      allocate(tt(np,np))
      do il=1,ll
        il3=3-il
        if(ll==1) il3=1
c --------- Real and imaginary parts separately ---------------
        do it=1,2
          if(it==1) tt=real(t_rs(:,:,il3))
          if(it==1) tt=imag(t_rs(:,:,il3))
          call dgemm('n','t',np,n*n,np,1.d0,tt,np,fif,n*n,0.d0,aa,np)
          call dgemm('n','n',n*n,n*n,np,1.d0,fif,n*n,aa,np,0.d0,bb,n*n)
          cc=0.d0
          do j1=1,n
            do i1=1,n
              do j=1,n
                do i=1,n
                  cc(i,j)=cc(i,j)+gtau(i1,j1,i_tau)*bb(j1,j,i1,i)
                enddo
              enddo
            enddo
          enddo
          if(it==1) vx(:,:,il)=vx(:,:,il)+cc
          if(it==2) vx(:,:,il)=vx(:,:,il)+(0.d0,1.d0)*cc
        enddo
      enddo   !! over il
      deallocate(tt,aa,bb,cc)
      end
