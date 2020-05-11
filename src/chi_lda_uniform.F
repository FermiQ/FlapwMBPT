      subroutine chi_lda_uniform(ch,xi0_nu,xi_nu,
     &                           n_ix0,ind_ix0,n_ch,eps1,eps10)
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
      real*8, intent(inout) :: xi0_nu(n_ch),xi_nu(n_ch),eps10,eps1
      integer :: n_pbt,i,info,j,n,nn
      real*8 :: pi2a
      complex*16 :: cc
      integer, allocatable :: ipiv(:)
      real*8, allocatable :: d(:),v_tmpr(:,:)
      complex*16, allocatable :: chi0(:,:),
     &                           work(:,:),thet_q(:,:,:),v_tmp(:,:),
     &                           chi(:,:),ixc0(:,:),tmp(:,:),tmp1(:,:),
     &                           a(:,:)
      pi2a=(pi+pi)/par
      n_pbt=n_pbmt+nplwgw(1)
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
      allocate(ipiv(nn))
      allocate(work(nn,nn))
      allocate(chi0(nn,nn))
      allocate(tmp(n_pbt,n_pbt))
      allocate(tmp1(n_pbmt,n_pbmt))
      allocate(chi(nn,nn))
      allocate(ixc0(nn,nn))
      allocate(a(nn,nn))
      allocate(d(nn))
c -------------------------------------------------------------------
      allocate(thet_q(n_pbt,n_pbt,n_ixc_0))
      call theta_q_ss(1,n_pbt,thet_q)
      if(ch=='___00'.or.ch=='00_ZZ') then
        allocate(v_tmp(n_pbt,n_pbt))
        v_tmp=(0.d0,0.d0)
        if(me_k==0) then
          if(ncmpl==1) then
            allocate(v_tmpr(n_pbt,n_pbt))
            call v_coul_full_r(nplwgw(1),v_tmpr,1,0)
            call pb_c_from_r(n_pbt,n_pbt,n_pbt,v_tmp,v_tmpr,pnt(1,1))
            deallocate(v_tmpr)
          else
            call v_coul_full(nplwgw(1),v_tmp,1,0)
          endif
        endif
        if(nproc_k/=1) call dgop(v_tmp,2*n_pbt**2,'  +',comm_k)
        thet_q(:,:,1)=thet_q(:,:,1)+v_tmp
      endif
      if(n==1.or.ch=='___XY') then
        j=ind_ix0(1)
        ixc0(:,:)=thet_q(:,:,j)
        if(ch=='___XY') ixc0(:,:)=ixc0(:,:)/2 
      else
        j=ind_ix0(1)
        ixc0(1:nn/2,1:nn/2)=thet_q(:,:,j)
        j=ind_ix0(2)
        ixc0(nn/2+1:nn,nn/2+1:nn)=thet_q(:,:,j)
        j=ind_ix0(3)
        ixc0(nn/2+1:nn,1:nn/2)=thet_q(:,:,j)
        ixc0(1:nn/2,nn/2+1:nn)=thet_q(:,:,j)
      endif
      deallocate(thet_q)
c --------------------------------------------------------------------
      if(nn==n_pbt) then
        if(ch=='___XY') then
          call chi0_lda_uniform(1,2,chi0,n_pbt)
          chi0=4.d0*chi0
        else if(ch=='___00'.or.ch=='___ZZ') then
          call chi0_lda_uniform(1,1,chi0,n_pbt)
          chi0=2.d0*chi0
        endif
      else if(nn>n_pbt) then
      endif
c  ------- Enhanced Susceptibility ------------------------------
      call zgemm('n','n',nn,nn,nn,(-1.d0,0.d0),chi0,nn,ixc0,nn,
     &           (0.d0,0.d0),work,nn)
      do i=1,nn
        work(i,i)=(1.d0,0.d0)+work(i,i)
      enddo
      chi=chi0
      call zgesv(nn,nn,work,nn,ipiv,chi,nn,info)
      if(me_k==0) then    
        call pw_from_pb(1,n_pbt,nn,chi0(1,1),cc,3)
        xi0_nu(1)=cc
        call pw_from_pb(1,n_pbt,nn,chi(1,1),cc,3)
        xi_nu(1)=cc
c ----------- Inverse Macroscopic Dielectric function ----------------
        if(ch=='___00'.or.ch=='00_ZZ') then
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &               n_pbt,chi0,nn,(0.d0,0.d0),work,nn)
	      do i=1,n_pbt
	        work(i,i)=(1.d0,0.d0)+work(i,i)
	      enddo
          call pw_from_pb(1,n_pbt,nn,work,cc,2)
          eps10=cc
          call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &               n_pbt,chi,nn,(0.d0,0.d0),work,nn)
	      do i=1,n_pbt
	        work(i,i)=(1.d0,0.d0)+work(i,i)
	      enddo
          call pw_from_pb(1,n_pbt,nn,work,cc,2)
          eps1=cc
        endif
        if(n_ch==3) then
          tmp=chi0(nn/2+1:nn,nn/2+1:nn)
          call pw_from_pb(1,n_pbt,n_pbt,tmp,cc,3)
          xi0_nu(2)=cc
          tmp=chi(nn/2+1:nn,nn/2+1:nn)
          call pw_from_pb(1,n_pbt,n_pbt,tmp,cc,3)
          xi_nu(2)=cc
          tmp=chi0(1:nn/2,nn/2+1:nn)
          call pw_from_pb(1,n_pbt,n_pbt,tmp,cc,3)
          xi0_nu(3)=cc
          tmp=chi(1:nn/2,nn/2+1:nn)
          call pw_from_pb(1,n_pbt,n_pbt,tmp,cc,3)
          xi_nu(3)=cc
        endif
      endif
      deallocate(a,d,ipiv,work,ixc0,chi0,tmp,tmp1,chi)
      if(ch=='___00'.or.ch=='00_ZZ') deallocate(v_tmp)
      end
