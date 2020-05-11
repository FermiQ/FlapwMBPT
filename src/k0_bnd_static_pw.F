      subroutine k0_bnd_static_pw(key,iq_line,iq,k0_pw,g_omega,g_x,ex)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: key,iq_line,iq
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c,nspin),
     &                      ex(nbndf_bnd,npnt_c,nspin)
      complex*16, intent(in) :: g_omega(nbndf_bnd,nbndf_bnd,0:n_omega,
     &                                  npnt_c,nspin)
	  complex*16, intent(out) :: k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,
     &                                 nqdiv_c,nspin)
      integer :: i_omega,k,kq,ispin,nn,n,i_tau,k0,nk,kq0,nkq,i,j,kf0,
     &           kfq0
      real*8 :: v(3),om,de
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:,:),kx(:,:),
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
	      v=pnt_c(:,k)-pnt_c(:,iq)
	      call zone1_number(v,rb0_c,ndiv_c,kq)
	      kq=index_k1_c(kq)     !! for K-Q
	      kq0=i_kref_c(kq)
          kfq0=k_a_from_c(kq0)
	      nkq=n_low_bnd(kfq0,ispin)
          allocate(tmp(nk,nkq))
          allocate(tmp1(nk,nkq,0:n_omega,2))
          allocate(kx(nk,nkq))
	      tmp1=(0.d0,0.d0)
	      do i_omega=0,n_omega
	        om=w_omega(i_omega)
c -------------- Negative omega (in=1) --------------------------------
	        call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),
     &	               psi_psi_pw(1,1,k,iq_line,ispin,ispin),n,
     &                 g_omega(1,1,i_omega,kq0,ispin),n,(0.d0,0.d0),
     &                 tmp,nk)
	        call zgemm('c','n',nk,nkq,nk,(-1.d0,0.d0),
     &	               g_omega(1,1,i_omega,k0,ispin),n,tmp,nk,
     &                 (0.d0,0.d0),tmp1(1,1,i_omega,1),nk)
c -------- Exchange part ----------------------------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-psi_psi_pw(i,j,k,iq_line,ispin,ispin)
     &                  /dcmplx(-ex(i,k0,ispin),-om)
     &                  /dcmplx(-ex(j,kq0,ispin),-om)
              enddo
            enddo
            tmp1(:,:,i_omega,1)=tmp1(:,:,i_omega,1)-kx
c -------------- Positive omega (in=2) --------------------------------
	        call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),
     &	               psi_psi_pw(1,1,k,iq_line,ispin,ispin),n,
     &                 g_omega(1,1,i_omega,kq0,ispin),n,(0.d0,0.d0),
     &                 tmp,nk)
	        call zgemm('n','n',nk,nkq,nk,(-1.d0,0.d0),
     &	               g_omega(1,1,i_omega,k0,ispin),n,tmp,nk,
     &                 (0.d0,0.d0),tmp1(1,1,i_omega,2),nk)
c -------- Exchange part -------------------------
            do j=1,nkq
              do i=1,nk
                kx(i,j)=-psi_psi_pw(i,j,k,iq_line,ispin,ispin)
     &                  /dcmplx(-ex(i,k0,ispin),om)
     &                  /dcmplx(-ex(j,kq0,ispin),om)
              enddo
            enddo
            tmp1(:,:,i_omega,2)=tmp1(:,:,i_omega,2)-kx
c --------- Transform to A and B form ---------------------------------
c --------- A=tmp1(+) + tmp1(-) ---------------------------------------
c --------- B=i{-tmp1(+) + tmp1(-)} -----------------------------------
            tmp=tmp1(:,:,i_omega,2)+tmp1(:,:,i_omega,1)
            tmp1(:,:,i_omega,1)=(0.d0,1.d0)*(-tmp1(:,:,i_omega,2)
     &                                       +tmp1(:,:,i_omega,1))
            tmp1(:,:,i_omega,2)=tmp
          enddo  !! over i_omega
c --------- Get A(tau) and B(tau) -------------------------------------
          allocate(a_tau(nk,nkq,0:n_tau/2))
          allocate(b_tau(nk,nkq,0:n_tau/2))
c ---------------------- A-form --------------------------------          
	      call dgemm('n','n',2*nk*nkq,n_tau/2+1,n_omega+1,1.d0,
     &	             tmp1(1,1,0,2),2*nk*nkq,tau_from_omega(0,0,1),
     &               n_omega+1,0.d0,a_tau,2*nk*nkq)
c ---------------------- B-form -------------------------------- 
	      call dgemm('n','n',2*nk*nkq,n_tau/2+1,n_omega+1,1.d0,
     &	             tmp1(1,1,0,1),2*nk*nkq,tau_from_omega(0,0,2),
     &               n_omega+1,0.d0,b_tau,2*nk*nkq)
c -------- Get tau and b-tau components of K0_PW ----------------------
          do i_tau=0,n_tau/2
            k0_pw(1:nk,1:nkq,i_tau,k,ispin)=
     &        0.25d0*(a_tau(:,:,i_tau)+b_tau(:,:,i_tau))
            k0_pw(1:nk,1:nkq,n_tau-i_tau,k,ispin)=
     &        0.25d0*(b_tau(:,:,i_tau)-a_tau(:,:,i_tau))
          enddo
          deallocate(a_tau,b_tau)
c --------- Add the asymptotic correction --------------------------
          do i_tau=0,n_tau
            do j=1,nkq
              do i=1,nk
                de=ex(j,kq0,ispin)-ex(i,k0,ispin)
                if(abs(de)<1.d-8) then
                  kx(i,j)=g_x(i,i_tau,k0,ispin)
     &              *(betta_t*g_x(i,n_tau,k0,ispin)+tau_mesh(i_tau))
                else
                  kx(i,j)=(g_x(i,i_tau,k0,ispin)
     &                    -g_x(j,i_tau,kq0,ispin))/de
                endif
                kx(i,j)=kx(i,j)
     &                 *psi_psi_pw(i,j,k,iq_line,ispin,ispin)
              enddo
            enddo
	        k0_pw(1:nk,1:nkq,i_tau,k,ispin)=
     &	          k0_pw(1:nk,1:nkq,i_tau,k,ispin)+kx
          enddo
          deallocate(tmp,tmp1,kx)
        enddo   !! over k
        if(key==1.and.nspin==2.and.ispin==1) then
          k0_pw(:,:,:,:,ispin)=-k0_pw(:,:,:,:,ispin)
        endif
      enddo   !! over ispin
      end
