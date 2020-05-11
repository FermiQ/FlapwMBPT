      subroutine k0_bnd(ind_nu,iq,k0_pw,g_om_nu,ps3,s2,g_x,ex,nom)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,iq,s2,nom
      real*8, intent(in) :: g_x(nbndf_bnd,0:n_tau,npnt_c,nspin),
     &                      ex(nbndf_bnd,npnt_c,nspin)
      complex*16, intent(in) :: g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c,nspin),
     &                          ps3(nbndf_bnd,nbndf_bnd,n_pbtot_red,
     &                              nqdiv_c,nspin)
	  complex*16, intent(out) :: k0_pw(nbndf_bnd,nbndf_bnd,0:n_tau,
     &                                 2,nqdiv_c,nspin)
      integer :: i_omega,k,kq,ispin,nn,in,n,i_tau,k0,nk,kq0,i_nu,
     &           nkq,i,j,kf0,kfq0
      real*8 :: v(3),om,omnu,de
      complex*16, allocatable :: tmp(:,:),tmp1(:,:,:),kx(:,:)
      i_nu=me_t*ndim3_nu+ind_nu-1
      n=nbndf_bnd
      nn=n*n
      allocate(tmp(n,n))
      allocate(tmp1(n,n,n_omega1_max))
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
          allocate(kx(nk,nkq))
	      do in=1,2
	        tmp1=(0.d0,0.d0)
	        do i_omega=1,nom
	          om=w_om_adapt_nu(i_omega,i_nu)
              omnu=om-w_nu(i_nu)
	          if(in==1) then
	            call zgemm('n','c',nk,nkq,nkq,(1.d0,0.d0),
     &                     ps3(1,1,s2,k,ispin),n,
     &                     g_om_nu(1,1,1,i_omega,kq0,ispin),n,
     &                     (0.d0,0.d0),tmp,n)
	            call zgemm('c','n',nk,nkq,nk,(1.d0,0.d0),
     &                     g_om_nu(1,1,2,i_omega,k0,ispin),n,
     &                     tmp,n,(0.d0,0.d0),tmp1(1,1,i_omega),n)
c -------- Exchange part -------------------------
                do j=1,nkq
                  do i=1,nk
                    kx(i,j)=-ps3(i,j,s2,k,ispin)
     &                     /dcmplx(-ex(i,k0,ispin),-omnu)
     &                     /dcmplx(-ex(j,kq0,ispin),-om)
                  enddo
                enddo
	          else if(in==2) then
	            call zgemm('n','n',nk,nkq,nkq,(1.d0,0.d0),
     &                     ps3(1,1,s2,k,ispin),n,
     &                     g_om_nu(1,1,2,i_omega,kq0,ispin),n,
     &                     (0.d0,0.d0),tmp,n)
	            call zgemm('n','n',nk,nkq,nk,(1.d0,0.d0),
     &                     g_om_nu(1,1,1,i_omega,k0,ispin),n,
     &                     tmp,n,(0.d0,0.d0),tmp1(1,1,i_omega),n)
c -------- Exchange part -------------------------
                do j=1,nkq
                  do i=1,nk
                    kx(i,j)=-ps3(i,j,s2,k,ispin)
     &                     /dcmplx(-ex(i,k0,ispin),om)
     &                     /dcmplx(-ex(j,kq0,ispin),omnu)
                  enddo
                enddo
              endif
              tmp1(:,:,i_omega)=-tmp1(:,:,i_omega)
c -------- Subtract Exchange part -------------------------
              tmp1(1:nk,1:nkq,i_omega)=tmp1(1:nk,1:nkq,i_omega)-kx
              if(in==1) tmp1(:,:,i_omega)=conjg(tmp1(:,:,i_omega))
            enddo  !! over i_omega
            call from_omega_nu_to_tau_nu_ab(nn,ind_nu,tmp1,
     &                                      k0_pw(1,1,0,in,k,ispin),in)
c --------- Add the asymptotic correction --------------------------
            do i_tau=0,n_tau
              do j=1,nkq
                do i=1,nk
                  de=ex(j,kq0,ispin)-ex(i,k0,ispin)
                  if(i_nu==0.and.abs(de)<1.d-8) then
                    kx(i,j)=0.5d0*g_x(i,i_tau,k0,ispin)
     &                     *(betta_t*g_x(i,n_tau,k0,ispin)
     &                              +tau_mesh(i_tau))
                  else if(in==1) then
                    kx(i,j)=g_x(i,i_tau,k0,ispin)
     &                     /dcmplx(de,w_nu(i_nu))
                  else if(in==2) then
                    kx(i,j)=-g_x(j,i_tau,kq0,ispin)
     &                      /dcmplx(de,w_nu(i_nu))
                  endif
                  kx(i,j)=kx(i,j)*ps3(i,j,s2,k,ispin)
                enddo
              enddo
              k0_pw(1:nk,1:nkq,i_tau,in,k,ispin)=
     &           k0_pw(1:nk,1:nkq,i_tau,in,k,ispin)+kx
            enddo
          enddo   !! over in
          deallocate(kx)
        enddo   !! over k
      enddo   !! over ispin
      deallocate(tmp,tmp1)
      end
