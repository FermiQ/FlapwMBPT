      subroutine sigma_g_lambda_w_02(n,np,nom,ind_nu,lambda_dyn,g_om,
     &                               w3_nu_d,aw,bw)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: n,np,nom,ind_nu
      real*8, intent(in) :: w3_nu_d(np,n,n)
      complex*16, intent(in) :: lambda_dyn(n,n,np,nom,2),
     &                          g_om(n,n,nom,2)
	complex*16, intent(out) :: aw(n,n,0:n_tau),bw(n,n,0:n_tau)
      integer :: i_omega,i_tau,ii,j,in,k
      complex*16, allocatable :: tmp(:,:,:,:),tmp1(:,:,:,:),a(:,:,:),
     &                           b(:,:)
      allocate(tmp(n,n,np,nom))
      allocate(tmp1(n,n,np,0:n_tau))
      allocate(a(n,np,n))
      allocate(b(n,n))
c ------------------------------------------------------------------
      aw=(0.d0,0.d0)
      bw=(0.d0,0.d0)
      do in=1,2
	  do i_omega=1,nom
	    do j=1,np
	      if(in==1) then
	        b=conjg(lambda_dyn(:,:,j,i_omega,1))
	        call zgemm('n','n',n,n,n,(1.d0,0.d0),g_om(1,1,i_omega,2),
     &	                 n,b,n,(0.d0,0.d0),tmp(1,1,j,i_omega),n)
	      else if(in==2) then
	        call zgemm('n','n',n,n,n,(1.d0,0.d0),g_om(1,1,i_omega,1),
     &	                 n,lambda_dyn(1,1,j,i_omega,2),n,(0.d0,0.d0),
     &	                 tmp(1,1,j,i_omega),n)
            endif
          enddo
        enddo  !! over i_omega
        call from_omega_nu_to_tau_nu_ab(n*n*np,ind_nu,tmp,tmp1,in)
        do i_tau=0,n_tau
          do j=1,n
            do k=1,n
              do ii=1,np
                a(j,ii,k)=tmp1(k,j,ii,i_tau)
              enddo
            enddo
          enddo
	    call dgemm('n','n',2*n,n,n*np,1.d0,a,2*n,w3_nu_d,np*n,0.d0,b,
     &	             2*n)
          if(in==1) aw(:,:,i_tau)=transpose(b)
          if(in==2) bw(:,:,i_tau)=transpose(b)
        enddo  !! over i_tau
      enddo  !! over in
      deallocate(tmp1,tmp,a,b)
      end