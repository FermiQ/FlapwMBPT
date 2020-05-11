      subroutine tsig_r_to_k_mi(t_rs,t_pw,k1,phf,ax,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: phf(nr_full_red,nqdiv_c),
     &                          t_rs(n_pbmt_red,nrel*nr_full_red,
     &                               nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c)
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      integer :: kq,kq0,nkq,iq,iq0,kqa0,kqa,irl,ind,ir,i,i0,ii
      real*8 :: v(3)
      complex*16 :: cc
      complex*16, allocatable :: tm(:,:),t(:,:),t1(:,:)
      allocate(tm(n_pbmt_red,nbndf_bnd))
      allocate(t(nr_full_red,n_pbmt_red))
      allocate(t1(n_pbmt_red,nrel*nplw_gw_red))
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        v=pnt_c(:,k1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K'-Q
        kqa=k_a_from_c(kq)
        kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
        nkq=n_low_bnd(kqa0,ispin)
        do irl=1,nrel
          ind=(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=conjg(phf(ir,kq))
            t(ir,:)=cc*t_rs(:,ind+ir,iq)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),n_pbmt_red,
     &              t,-1)
          ii=(irl-1)*nplwgw_red(kqa0)
          do i=1,nplwgw_red(kqa0)
            i0=indpw_gw_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(:,ii+i)=t(ind,:)
          enddo
        enddo   !! over irl
        t1=t1/nr_full_red*sqrt(amega)
        t1=conjg(t1)
        call zgemm('n','n',n_pbmt_red,nkq,nrel*nplwgw_red(kqa0),
     &             (1.d0,0.d0),t1,n_pbmt_red,ax(1,1,kq),
     &             nrel*nplw_gw_red,(0.d0,0.d0),tm,n_pbmt_red)
        tm=conjg(tm)
        t_pw(1:n_pbmt_red,1:nkq,iq)=t_pw(1:n_pbmt_red,1:nkq,iq)
     &                             +tm(:,1:nkq)
      enddo   !! over iq
      deallocate(tm,t1,t)
      end
      
      
      subroutine tsig_r_to_k_im(t_rs,t_pw,k1,z_red,phf,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: phf(nr_full_red,nqdiv_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          t_rs(nr_full_red,nfun_red,nqdiv_c)
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      integer :: kq,kq0,nkq,kqa0,iq,iq0,iqa0,iqa,np,ir,i,i0,ind,j,
     &           jnd,id,gbs_number
      real*8 :: v(3)
      complex*16 :: cc
      real*8, allocatable :: sr(:,:),tr(:,:),tr1(:,:)
      complex*16, allocatable :: tm(:,:),tt(:,:),t(:,:),t1(:,:),s(:,:)
      allocate(tm(nfun_red,nbndf_bnd))
      allocate(t(nr_full_red,nfun_red))
      allocate(t1(nplw_gw_red,nfun_red))
      allocate(tt(nplw_gw_red,nbndf_bnd))
      if(complex_ro) then
        allocate(s(nplw_gw_red,nplw_gw_red))
      else
        allocate(sr(nplw_gw_red,nplw_gw_red))
        allocate(tr(nplw_gw_red,nbndf_bnd))
        allocate(tr1(nplw_gw_red,nbndf_bnd))
      endif
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
        iqa=k_a_from_c(iq)
        v=pnt_c(:,k1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K'-Q
        kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
        nkq=n_low_bnd(kqa0,ispin)
        np=nplwgw_red(iqa0)
        do ir=1,nr_full_red
          cc=conjg(phf(ir,iq))
          t(ir,:)=cc*t_rs(ir,:,iq)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nfun_red,t,-1)
        do i=1,np
          i0=indpw_gw_red(i,iqa)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t1(i,:)=t(ind,:)
        enddo
        t1=t1/nr_full_red*sqrt(amega)
        tm(:,1:nkq)=conjg(z_red(:,1:nkq,kq))
        call zgemm('n','n',np,nkq,nfun_red,(1.d0,0.d0),t1,nplw_gw_red,
     &             tm,nfun_red,(0.d0,0.d0),tt,nplw_gw_red)
        if(complex_ro) then
          do j=1,np
            jnd=indpw_gw_red(j,iqa)
            do i=1,np
              ind=indpw_gw_red(i,iqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              s(i,j)=dcmplx(sovr(id),sovi(id))
            enddo
          enddo
          call zgemm('n','n',np,nkq,np,(1.d0,0.d0),s,nplw_gw_red,tt,
     &               nplw_gw_red,(1.d0,0.d0),
     &               t_pw(n_pbmt_red+1,1,iq),n_pbtot_red)
        else
          do j=1,np
            jnd=indpw_gw_red(j,iqa)
            do i=1,np
              ind=indpw_gw_red(i,iqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              sr(i,j)=sovr(id)
            enddo
          enddo
          tr(1:np,1:nkq)=real(tt(1:np,1:nkq))
          call dgemm('n','n',np,nkq,np,1.d0,sr,nplw_gw_red,tr,
     &               nplw_gw_red,0.d0,tr1,nplw_gw_red)
          do j=1,nkq
            do i=1,np
              t_pw(n_pbmt_red+i,j,iq)=t_pw(n_pbmt_red+i,j,iq)+tr1(i,j)
            enddo
          enddo
          tr(1:np,1:nkq)=imag(tt(1:np,1:nkq))
          call dgemm('n','n',np,nkq,np,1.d0,sr,nplw_gw_red,tr,
     &               nplw_gw_red,0.d0,tr1,nplw_gw_red)
          do j=1,nkq
            do i=1,np
              t_pw(n_pbmt_red+i,j,iq)=t_pw(n_pbmt_red+i,j,iq)
     &                               +(0.d0,1.d0)*tr1(i,j)
            enddo
          enddo
        endif
      enddo   !! over iq
      deallocate(tm,t,t1,tt)
      if(complex_ro) deallocate(s)
      if(.not.complex_ro) deallocate(sr,tr,tr1)
      end
      
      
      
      
      subroutine tsig_r_to_k_ii(t_rs,t_pw,k1,ax,phf,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: phf(nr_full_red,nqdiv_c),
     &                          ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c),
     &                          t_rs(nr_full_red,nrel*nr_full_red,
     &                               nqdiv_c)
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      integer :: kq,kq0,nkq,iq,iq0,iqa0,iqa,kqa0,kqa,np,np1,
     &           irl,ind,ir,ii,i,i0,gbs_number,id,j,jnd
      real*8 :: v(3)
      complex*16 :: cc
      real*8, allocatable :: sr(:,:),tr(:,:),tr1(:,:)
      complex*16, allocatable :: tm(:,:),s(:,:),t(:,:),t1(:,:),t2(:,:)
      allocate(tm(nplw_gw_red,nbndf_bnd))
      allocate(t(nr_full_red,nr_full_red))
      allocate(t1(nr_full_red,nrel*nplw_gw_red))
      allocate(t2(nplw_gw_red,nrel*nplw_gw_red))
      if(complex_ro) then
        allocate(s(nplw_gw_red,nplw_gw_red))
      else
        allocate(sr(nplw_gw_red,nplw_gw_red))
        allocate(tr(nplw_gw_red,nbndf_bnd))
        allocate(tr1(nplw_gw_red,nbndf_bnd))
      endif
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
        iqa=k_a_from_c(iq)
        v=pnt_c(:,k1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,kq)
        kq=index_k1_c(kq)     !! for K'-Q
        kq0=i_kref_c(kq)
        kqa0=k_a_from_c(kq0)
        kqa=k_a_from_c(kq)
        nkq=n_low_bnd(kqa0,ispin)
        np=nplwgw_red(iqa0)
        np1=nplwgw_red(kqa0)
        do irl=1,nrel
          ind=(irl-1)*nr_full_red
          do ir=1,nr_full_red
            cc=conjg(phf(ir,kq))
            t(ir,:)=cc*t_rs(:,ind+ir,iq)
          enddo
          call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nr_full_red,
     &              t,-1)
          ii=(irl-1)*np1
          do i=1,np1
            i0=indpw_gw_red(i,kqa)
            call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
            t1(:,ii+i)=t(ind,:)
          enddo
        enddo   !! over irl
        t1=t1/nr_full_red*sqrt(amega)
        do ir=1,nr_full_red
          cc=conjg(phf(ir,iq))
          t1(ir,:)=cc*t1(ir,:)
        enddo
        call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nrel*np1,t1,-1)
        do i=1,np
          i0=indpw_gw_red(i,iqa)
          call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
          t2(i,:)=t1(ind,:)
        enddo
        t2=t2/nr_full_red*sqrt(amega)
        t2=conjg(t2)
        call zgemm('n','n',np,nkq,nrel*np1,(1.d0,0.d0),t2,nplw_gw_red,
     &             ax(1,1,kq),nrel*nplw_gw_red,(0.d0,0.d0),tm,
     &             nplw_gw_red)
        tm=conjg(tm)
        if(complex_ro) then
          do j=1,np
            jnd=indpw_gw_red(j,iqa)
            do i=1,np
              ind=indpw_gw_red(i,iqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              s(i,j)=dcmplx(sovr(id),sovi(id))
            enddo
          enddo
          call zgemm('n','n',np,nkq,np,(1.d0,0.d0),s,nplw_gw_red,tm,
     &               nplw_gw_red,(1.d0,0.d0),
     &               t_pw(n_pbmt_red+1,1,iq),n_pbtot_red)
        else
          do j=1,np
            jnd=indpw_gw_red(j,iqa)
            do i=1,np
              ind=indpw_gw_red(i,iqa)
              v=gbs(:,jnd)-gbs(:,ind)
              id=gbs_number(v)
              sr(i,j)=sovr(id)
            enddo
          enddo
          tr(1:np,1:nkq)=real(tm(1:np,1:nkq))
          call dgemm('n','n',np,nkq,np,1.d0,sr,nplw_gw_red,tr,
     &               nplw_gw_red,0.d0,tr1,nplw_gw_red)
          do j=1,nkq
            do i=1,np
              t_pw(n_pbmt_red+i,j,iq)=t_pw(n_pbmt_red+i,j,iq)+tr1(i,j)
            enddo
          enddo
          tr(1:np,1:nkq)=imag(tm(1:np,1:nkq))
          call dgemm('n','n',np,nkq,np,1.d0,sr,nplw_gw_red,tr,
     &               nplw_gw_red,0.d0,tr1,nplw_gw_red)
          do j=1,nkq
            do i=1,np
              t_pw(n_pbmt_red+i,j,iq)=t_pw(n_pbmt_red+i,j,iq)
     &                               +(0.d0,1.d0)*tr1(i,j)
            enddo
          enddo
        endif
      enddo   !! over iq
      deallocate(tm,t,t1,t2)
      if(complex_ro) deallocate(s)
      if(.not.complex_ro) deallocate(sr,tr,tr1)
      end
