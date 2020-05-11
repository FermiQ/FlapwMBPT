      subroutine k_t_sig_sing(ispin,k1,k_pw,k_sing,gctau,t_pw,t_sing,
     &                        ps3)
      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ispin,k1
      complex*16, intent(in) :: gctau(nbndf_bnd,nbndf_bnd,npnt_c),
     &                          k_pw(nbndf_bnd,n_pbtot_red,nqdiv_c),
     &                          ps3(n_pbtot_red,nbndf_bnd,nbndf_bnd,
     &                              nqdiv_c),
     &                          k_sing(nbndf_bnd)
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      complex*16, intent(out) :: t_sing(nbndf_bnd)
      integer :: kf,n,iq,iq0,iqa0,kq,kq0,kqf,n1,np,i,j,ii,k,k0,ka
      real*8 :: w,v(3)
      complex*16, allocatable :: t(:,:),t1(:,:)
      w=-1.d0/sqrt(amega)/nqdiv_c
c ----- Contribution to regular T ----------------------------------
      kf=k_a_from_c(k1)
      n=n_low_bnd(kf,ispin)
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
        v=pnt_c(:,k1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)    !! for K-Q
        kq0=i_kref_c(kq)
        kqf=k_a_from_c(kq0)
        n1=n_low_bnd(kqf,ispin)
        np=n_pbmt_red+nplwgw_red(iqa0)
        allocate(t(np,n1))
        allocate(t1(np,n1))
        t=(0.d0,0.d0)
        do j=1,n
          do i=1,n1
            do ii=1,np
              t(ii,i)=t(ii,i)+conjg(ps3(ii,i,j,kq))*k_sing(j)
            enddo
          enddo
        enddo
        call zgemm('n','t',np,n1,n1,(1.d0,0.d0),t,np,gctau(1,1,kq0),
     &             nbndf_bnd,(0.d0,0.d0),t1,np)
        do j=1,n1
          do i=1,np
            t_pw(i,j,iq)=t_pw(i,j,iq)+t1(i,j)*w
          enddo
        enddo
        deallocate(t,t1)
      enddo
c ----- Contribution to singular T ---------------------------------
      t_sing=0.d0
      do k=1,nqdiv_c
        k0=i_kref_c(k)
        ka=k_a_from_c(k0)
        v=pnt_c(:,k1)-pnt_c(:,k)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)    !! for K-K'
        kq0=i_kref_c(kq)
        kqf=k_a_from_c(kq0)
        n1=n_low_bnd(ka,ispin)
        np=n_pbmt_red+nplwgw_red(kqf)
        allocate(t(n1,n1))
        do ii=1,n
          call zgemm('n','n',n1,n1,np,(1.d0,0.d0),k_pw(1,1,k),nbndf_bnd,
     &               ps3(1,1,ii,k),n_pbtot_red,(0.d0,0.d0),t,n1)
          do j=1,n1
            do i=1,n1
              t_sing(ii)=t_sing(ii)+w*gctau(i,j,k0)*t(j,i)
            enddo
          enddo
        enddo
        deallocate(t)
      enddo
      w=w/sqrt(amega)
      do ii=1,n
        do i=1,n
          t_sing(ii)=t_sing(ii)+w*gctau(ii,i,k1)*k_sing(i)
        enddo
      enddo
      end
