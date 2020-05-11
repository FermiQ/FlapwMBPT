      subroutine chi_lda_static(ch,xi0_nu,xi_nu,n_ix0,ind_ix0,n_ch,eps1,
     &                          eps10)
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
      real*8, intent(inout) :: xi0_nu(n_q_long,n_ch),
     &                         xi_nu(n_q_long,n_ch),eps10(n_q_long),
     &                         eps1(n_q_long)
      integer :: n_pbt,i,n,nn,nq0,iq,ind_q,iq0,ind,ii
      real*8 :: pi2a
      complex*16 :: cc
      real*8, allocatable :: v_tmpr(:,:)
      complex*16, allocatable :: chi0(:,:),chi(:,:),v_tmp(:,:),
     &                           thet_q(:,:,:),work(:,:),tmp(:,:)
      pi2a=(pi+pi)/par
      do nq0=1,n_q0_long
        iq=q0_long_in_npnt(nq0)
        n_pbt=n_pbmt+nplwgw(iq)
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
        call theta_q_ss(iq,n_pbt,thet_q)
        if(ch=='___00'.or.ch=='00_ZZ') then
          allocate(v_tmp(n_pbt,n_pbt))
          v_tmp=(0.d0,0.d0)
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
          if(nproc_k/=1) call dgop(v_tmp,2*n_pbt**2,'  +',comm_k)
          thet_q(:,:,1)=thet_q(:,:,1)+v_tmp
        endif
        allocate(chi0(nn,nn))
        allocate(chi(nn,nn))
        call chi_lda_static_1q(iq,n,nn,n_pbt,ch,n_ix0,ind_ix0,thet_q,
     &                         chi0,chi)
        deallocate(thet_q)
        allocate(tmp(n_pbt,n_pbt))
c -------------------------------------------------------------        
        do ii=1,size_q_long(nq0)
          ind=index_q_long(2,ii,nq0)  !! in Q_long list
          call pw_from_pb_long(nn,n_pbt,chi0,cc,ind,3)
          xi0_nu(ind,1)=cc
          call pw_from_pb_long(nn,n_pbt,chi,cc,ind,3)
          xi_nu(ind,1)=cc
c ----------- Inverse Macroscopic Dielectric function ----------------
          if(ch=='___00'.or.ch=='00_ZZ') then
            allocate(work(nn,nn))
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                 n_pbt,chi0,nn,(0.d0,0.d0),work,nn)
	        do i=1,n_pbt
	          work(i,i)=(1.d0,0.d0)+work(i,i)
	        enddo
            call pw_from_pb_long(nn,n_pbt,work,cc,ind,2)  
            eps10(ind)=cc
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,(1.d0,0.d0),v_tmp,
     &                 n_pbt,chi,nn,(0.d0,0.d0),work,nn)
	        do i=1,n_pbt
	          work(i,i)=(1.d0,0.d0)+work(i,i)
	        enddo
            call pw_from_pb_long(nn,n_pbt,work,cc,ind,2)  
            eps1(ind)=cc
            deallocate(work)
          endif
          if(n_ch==3) then
            tmp=chi0(nn/2+1:nn,nn/2+1:nn)
            call pw_from_pb_long(n_pbt,n_pbt,tmp,cc,ind,3)
            xi0_nu(ind,2)=cc
            tmp=chi(nn/2+1:nn,nn/2+1:nn)
            call pw_from_pb_long(n_pbt,n_pbt,tmp,cc,ind,3)
            xi_nu(ind,2)=cc
            tmp=chi0(1:nn/2,nn/2+1:nn)
            call pw_from_pb_long(n_pbt,n_pbt,tmp,cc,ind,3)
            xi0_nu(ind,3)=cc
            tmp=chi(1:nn/2,nn/2+1:nn)
            call pw_from_pb_long(n_pbt,n_pbt,tmp,cc,ind,3)
            xi_nu(ind,3)=cc
          endif
        enddo
        deallocate(chi0,chi,tmp)
        if(ch=='___00'.or.ch=='00_ZZ') deallocate(v_tmp)
      enddo   !! over nq0
      end
