      subroutine k_lin_vec
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i,j,k,ind
      real*8 :: pi2,pi2a,qb(3,3),qq0(3)
      pi2=pi+pi
      pi2a=pi2/par
      nn_k_int=n_k_int(1)*n_k_int(2)*n_k_int(3)
      allocate(k_int_vec(3,nn_k_int))
      do i=1,3
        qb(:,i)=qb0(:,i)/n_k_int(i)
      enddo
      qq0=0.5d0*(qb0(:,1)+qb0(:,2)+qb0(:,3))
      ind=0
      do k=1,n_k_int(3)
        do j=1,n_k_int(2)
          do i=1,n_k_int(1)
            ind=ind+1
            k_int_vec(:,ind)=pi2a*(qb(:,1)*(i-0.5d0)+qb(:,2)*(j-0.5d0)
     &                            +qb(:,3)*(k-0.5d0)-qq0)
          enddo
        enddo
      enddo
      end
