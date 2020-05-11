      subroutine chi_lda_dynamic(ch,wrx,xi0,xi,n_ix0,ind_ix0,n_ch,
     &                           eps1,eps10)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ix0,ind_ix0(n_ix0),n_ch
      real*8, intent(in) :: wrx(0:nrax_chi)
      complex*16, intent(inout) :: xi0(0:nrax_chi,n_ch),
     &                             xi(0:nrax_chi,n_ch),
     &                             eps10(0:nrax_chi),
     &                             eps1(0:nrax_chi)
      integer :: n_pbt,i,n,nn,iq,ind_i,ind_k,k,ii,info,i0,j,ind_q,iq0
      real*8 :: pi2a,v(3),v2
      integer, allocatable :: nd(:),ns(:),ipiv(:)
      real*8, allocatable :: vcou(:),v_tmpr(:,:)
      complex*16, allocatable :: chi0(:,:,:),chi(:,:,:),v_tmp(:,:),
     &                           thet_q(:,:,:),work(:,:),tmp(:,:),pw(:)
      pi2a=(pi+pi)/par
      iq=q_suscept
      n_pbt=nplwgw(iq)
      if(chi_basis=='PB') n_pbt=n_pbmt+nplwgw(iq)
      if(ch=='___00'.or.ch=='___ZZ') then
        n=1
        nn=n_pbt
      else if(ch=='00_ZZ') then
        n=2
        nn=2*n_pbt
      else if(ch=='___XY') then
        n=2
        nn=n_pbt
      endif
      allocate(thet_q(n_pbt,n_pbt,n_ixc_0))
      if(chi_basis=='PB') call theta_q_ss(iq,n_pbt,thet_q)
      if(chi_basis=='PW') call theta_q_pw(iq,n_pbt,thet_q)
      if(ch=='___00'.or.ch=='00_ZZ') then
        if(chi_basis=='PB') then
          allocate(v_tmp(n_pbt,n_pbt))
          do ind_q=1,ndim3_k(me_k+1)
            iq0=n3_mpi_k(me_k+1)+ind_q
            if(iq0==iq) then
              if(ncmpl==1) then
                allocate(v_tmpr(n_pbt,n_pbt))
                call v_coul_full_r(nplwgw(iq),v_tmpr,ind_q,0)
                call pb_c_from_r(n_pbt,n_pbt,n_pbt,v_tmp,v_tmpr,
     &                           pnt(1,iq0))
                deallocate(v_tmpr)
              else
                call v_coul_full(nplwgw(iq),v_tmp,ind_q,0)
              endif
            endif
          enddo
          thet_q(:,:,1)=thet_q(:,:,1)+v_tmp
        else if(chi_basis=='PW') then
          allocate(vcou(n_pbt))
          do i=1,n_pbt
            i0=indpw_gw(i,iq)
            v=pi2a*(pnt(:,iq)+gbs(:,i0))
            v2=dot_product(v,v)
            vcou(i)=8.d0*pi/v2
            thet_q(i,i,1)=thet_q(i,i,1)+vcou(i)
          enddo
        endif
      endif
      allocate(nd(nproc_t))
      allocate(ns(nproc_t))
      call size_shift_par(nrax_chi+1,nproc_t,nd,ns)
      allocate(chi0(nn,nn,nd(me_t+1)))
      allocate(chi(nn,nn,nd(me_t+1)))
      call chi_lda_dynamic_1q(iq,n,nn,n_pbt,ch,n_ix0,ind_ix0,thet_q,
     &                        chi0,chi,wrx,nd(me_t+1),ns(me_t+1))
      deallocate(thet_q)
c  ------- Bare Coulomb enhancement for Charge case ----------------
      if(ch=='___00') then
        allocate(ipiv(nn))
        allocate(work(nn,nn))
        do ind_i=1,nd(me_t+1)
          i=ns(me_t+1)+ind_i-1
          if(chi_basis=='PB') then
            call zgemm('n','n',nn,nn,nn,(-1.d0,0.d0),chi0(1,1,ind_i),nn,
     &                 v_tmp,nn,(0.d0,0.d0),work,nn)
          else if(chi_basis=='PW') then
            do j=1,nn
              work(:,j)=-chi0(:,j,ind_i)*vcou(j)
            enddo
          endif
          do i=1,nn
            work(i,i)=(1.d0,0.d0)+work(i,i)
          enddo
          call zgesv(nn,nn,work,nn,ipiv,chi0(1,1,ind_i),nn,info)
        enddo
        deallocate(work,ipiv)
      endif
      allocate(tmp(n_pbt,n_pbt))
c -------------------------------------------------------------
      xi0=(0.d0,0.d0)
      xi=(0.d0,0.d0)
      if(ch=='___00'.or.ch=='00_ZZ') then
        eps1=(0.d0,0.d0)
        eps10=(0.d0,0.d0)
      endif
      if(chi_basis=='PB') then
        allocate(pw(n_pbt))
        pw=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          if(k/=iq) cycle
          pw=pw_pb(1:n_pbt,ind_k)
        enddo
        if(nproc_k/=1) call dgop(pw,2*n_pbt,'  +',comm_k)
        do ind_i=1,nd(me_t+1)
          i=ns(me_t+1)+ind_i-1
          call pw_from_pb_1q(n_pbt,chi0(1,1,ind_i),xi0(i,1),3,pw)
          call pw_from_pb_1q(n_pbt,chi(1,1,ind_i),xi(i,1),3,pw)
c ----------- Inverse Macroscopic Dielectric function --------------
          if(ch=='___00'.or.ch=='00_ZZ') then
            allocate(work(nn,nn))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                 n_pbt,chi0(1,1,ind_i),nn,(0.d0,0.d0),work,nn)
            do ii=1,n_pbt
              work(ii,ii)=(1.d0,0.d0)+work(ii,ii)
            enddo
            call pw_from_pb_1q(n_pbt,work,eps10(i),2,pw)
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                 n_pbt,chi(1,1,ind_i),nn,(0.d0,0.d0),work,nn)
            do ii=1,n_pbt
              work(ii,ii)=(1.d0,0.d0)+work(ii,ii)
            enddo
            call pw_from_pb_1q(n_pbt,work,eps1(i),2,pw)
            deallocate(work)
          endif
          if(n_ch==3) then
            tmp=chi0(nn/2+1:nn,nn/2+1:nn,ind_i)
            call pw_from_pb_1q(n_pbt,tmp,xi0(i,2),3,pw)
            tmp=chi(nn/2+1:nn,nn/2+1:nn,ind_i)
            call pw_from_pb_1q(n_pbt,tmp,xi(i,2),3,pw)
            tmp=chi0(1:nn/2,nn/2+1:nn,ind_i)
            call pw_from_pb_1q(n_pbt,tmp,xi0(i,3),3,pw)
            tmp=chi(1:nn/2,nn/2+1:nn,ind_i)
            call pw_from_pb_1q(n_pbt,tmp,xi(i,3),3,pw)
          endif
        enddo
      else if(chi_basis=='PW') then
        do ind_i=1,nd(me_t+1)
          i=ns(me_t+1)+ind_i-1
          xi0(i,1)=chi0(1,1,ind_i)
          xi(i,1)=chi(1,1,ind_i)

c ----------- Inverse Macroscopic Dielectric function --------------
          if(ch=='___00'.or.ch=='00_ZZ') then
            allocate(work(nn,nn))
            do ii=1,n_pbt
              work(ii,:)=vcou(ii)*chi0(ii,:,ind_i)
            enddo
            do ii=1,n_pbt
              work(ii,ii)=(1.d0,0.d0)+work(ii,ii)
            enddo
            eps10(i)=work(1,1)
            do ii=1,n_pbt
              work(ii,:)=vcou(ii)*chi(ii,:,ind_i)
            enddo
            do ii=1,n_pbt
              work(ii,ii)=(1.d0,0.d0)+work(ii,ii)
            enddo
            eps1(i)=work(1,1)
            deallocate(work)
          endif
          if(n_ch==3) then
            xi0(i,2)=chi0(nn/2+1,nn/2+1,ind_i)
            xi(i,2)=chi(nn/2+1,nn/2+1,ind_i)
            xi0(i,3)=chi0(1,nn/2+1,ind_i)
            xi(i,3)=chi(1,nn/2+1,ind_i)
          endif
        enddo
      endif
      deallocate(chi0,chi,tmp)
      if(chi_basis=='PB') deallocate(pw)
      if(ch=='___00'.or.ch=='00_ZZ') then
        if(chi_basis=='PB') deallocate(v_tmp)
        if(chi_basis=='PW') deallocate(vcou)
      endif
      deallocate(nd,ns)
      if(nproc_t/=1) then
        call dgop(xi0,2*n_ch*(nrax_chi+1),'  +',comm_t)
        call dgop(xi,2*n_ch*(nrax_chi+1),'  +',comm_t)
        if(ch=='___00'.or.ch=='00_ZZ') then
          call dgop(eps1,2*(nrax_chi+1),'  +',comm_t)
          call dgop(eps10,2*(nrax_chi+1),'  +',comm_t)
        endif
      endif
      end
