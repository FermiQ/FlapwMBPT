      subroutine k0_bnd_uniform(key,k0_pw,g_omega,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c,nspin),
     &                      ex(nbndf_bnd,npnt_c,nspin)
      complex*16, intent(in) :: g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,
     &                                  npnt_c,nspin)
	  complex*16, intent(out) :: k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,
     &                                 nqdiv_c,nspin)
      integer :: i_omega,k,ispin,nn,n,i_tau,k0,nk,i,kf0
      real*8 :: om
      complex*16 :: kx
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:,:),
     &                           a_tau(:,:,:),b_tau(:,:,:)
      n=nbndf_bnd
      nn=n*n
c ------------------------------------------------------------------
      k0_pw=(0.d0,0.d0)
      do ispin=1,nspin
        do k=1,nqdiv_c
	      k0=i_kref_c(k)
          kf0=k_a_from_c(k0)
	      nk=n_low_bnd(kf0,ispin)
          allocate(tmp(nk,nk))
          allocate(tmp1(nk,nk,0:n_omega,2))
	      tmp1=(0.d0,0.d0)
	      do i_omega=0,n_omega
	        om=w_omega(i_omega)
	        call zgemm('n','n',nk,nk,nk,(-1.d0,0.d0),
     &	               g_omega(1,1,i_omega,k0,ispin),n,
     &                 g_omega(1,1,i_omega,k0,ispin),n,(0.d0,0.d0),
     &                 tmp1(1,1,i_omega,2),nk)
c -------- Exchange part -------------------------
            do i=1,nk
              kx=-(1.d0,0.d0)/dcmplx(-ex(i,k0,ispin),om)**2
              tmp1(i,i,i_omega,2)=tmp1(i,i,i_omega,2)-kx
            enddo
            tmp1(:,:,i_omega,1)=conjg(transpose(tmp1(:,:,i_omega,2)))
c --------- Transform to A and B form ---------------------------------
c --------- A=tmp1(+) + tmp1(-) ---------------------------------------
c --------- B=i{-tmp1(+) + tmp1(-)} -----------------------------------
            tmp=tmp1(:,:,i_omega,2)+tmp1(:,:,i_omega,1)
            tmp1(:,:,i_omega,1)=(0.d0,1.d0)*(-tmp1(:,:,i_omega,2)
     &                                       +tmp1(:,:,i_omega,1))
            tmp1(:,:,i_omega,2)=tmp
          enddo  !! over i_omega
c --------- Get A(tau) and B(tau) -------------------------------------
          allocate(a_tau(nk,nk,0:n_tau/2))
          allocate(b_tau(nk,nk,0:n_tau/2))
c ---------------------- A-form --------------------------------          
	      call dgemm('n','n',2*nk*nk,n_tau/2+1,n_omega+1,1.d0,
     &	             tmp1(1,1,0,2),2*nk*nk,tau_from_omega(0,0,1),
     &               n_omega+1,0.d0,a_tau,2*nk*nk)
c ---------------------- B-form -------------------------------- 
	      call dgemm('n','n',2*nk*nk,n_tau/2+1,n_omega+1,1.d0,
     &	             tmp1(1,1,0,1),2*nk*nk,tau_from_omega(0,0,2),
     &               n_omega+1,0.d0,b_tau,2*nk*nk)
c -------- Get tau and b-tau components of K0_PW ----------------------
          do i_tau=0,n_tau/2
            k0_pw(1:nk,1:nk,i_tau,k,ispin)=
     &        0.25d0*(a_tau(:,:,i_tau)+b_tau(:,:,i_tau))
            k0_pw(1:nk,1:nk,n_tau-i_tau,k,ispin)=
     &        0.25d0*(b_tau(:,:,i_tau)-a_tau(:,:,i_tau))
          enddo
          deallocate(a_tau,b_tau)
c --------- Add the asymptotic correction --------------------------
          do i_tau=0,n_tau
            do i=1,nk
              kx=g_x(i,i_tau,k0,ispin)
     &          *(betta_t*g_x(i,n_tau,k0,ispin)+tau_mesh(i_tau))
	          k0_pw(i,i,i_tau,k,ispin)=k0_pw(i,i,i_tau,k,ispin)+kx
            enddo
          enddo
          deallocate(tmp,tmp1)
        enddo   !! over k
        if(key==1.and.nspin==2.and.ispin==1) then
          k0_pw(:,:,:,:,ispin)=-k0_pw(:,:,:,:,ispin)
        endif
      enddo   !! over ispin
      k0_pw=k0_pw/sqrt(amega)
      end
