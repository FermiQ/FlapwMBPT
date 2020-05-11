      subroutine t_matr_inf(v4_1,a,ii1,t)
c     ii1 - the order (in W) the summation starts with      
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ii1
      complex*16, intent(in) :: a(n_lem2_red,n_lem2_red,ndim3_nu,
     &                            ndimc_kk(me_k+1)),
     &                          v4_1(maxel_red**2,maxel_red**2,natom)
      complex*16, intent(out) :: t(n_lem2_red,n_lem2_red,ndim3_nu,
     &                             ndimc_kk(me_k+1))
      integer :: mx2,iatom,isort,ind2,n,nn,ind_nu,ind_q,j,i,info
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: vc(:,:)
      mx2=maxel_red**2
      allocate(vc(n_lem2_red,n_lem2_red))
      allocate(ipiv(n_lem2_red))
      t=(0.d0,0.d0)
      do ind_nu=1,ndim3_nu
        do ind_q=1,ndimc_kk(me_k+1)
c --- Right part: (WL)^ii1*W ----------------------  
	      do iatom=1,natom
	        isort=is(iatom)
	        n=lfunm_red(isort)
	        nn=n*n
            ind2=io_lem2_red(iatom)-1
            do j=1,nn
              do i=1,nn
                t(ind2+i,ind2+j,ind_nu,ind_q)=v4_1(i,j,iatom)
              enddo
            enddo
          enddo
          do i=2,ii1
	        call zgemm('n','n',n_lem2_red,n_lem2_red,n_lem2_red,
     &	               (1.d0,0.d0),a(1,1,ind_nu,ind_q),n_lem2_red,
     &                 t(1,1,ind_nu,ind_q),n_lem2_red,(0.d0,0.d0),vc,
     &                 n_lem2_red)
            t(:,:,ind_nu,ind_q)=vc
          enddo
c --- Left part: 1 - WL ----------------------  
          vc=-a(:,:,ind_nu,ind_q)
          do i=1,n_lem2_red
            vc(i,i)=(1.d0,0.d0)+vc(i,i)
          enddo
          call zgesv(n_lem2_red,n_lem2_red,vc,n_lem2_red,ipiv,
     &               t(1,1,ind_nu,ind_q),n_lem2_red,info)
        enddo
      enddo
      deallocate(vc,ipiv)
      end
