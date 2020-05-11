      subroutine p_ii_q_mem_4r(pg1,nd_b,n0_b,ind_tau,p)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nd_b,n0_b,ind_tau
      real*8, intent(in) :: pg1(nplw_gw,nd_b,ndim3_k(me_k+1))
      real*8, intent(inout) :: p(nd_b_pbint(me_b+1),ndim3_tn,
     &                           ndim3_k(me_k+1))
      integer :: k,npw,i,ind_k,j,ij,n2,n1,n3
      real*8, allocatable :: tmp(:,:),tmp1(:,:),s(:,:)
      n1=n0_b_pbint(me_b+1)
      n3=n1+nd_b_pbint(me_b+1)
      allocate(s(nplw_gw,nplw_gw))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        npw=nplwgw(k)
c ------ We form the overlap matrix --------------------------------
        call s_overlap_r(0,s,nplw_gw,npw,indpw_gw(1,k))
        n2=min(nd_b,npw-n0_b)
        allocate(tmp1(npw,npw))
        tmp1=0.d0
        if(n2>0) then
          allocate(tmp(npw,n2))
          call dgemm('t','n',npw,n2,npw,1.d0,s,nplw_gw,pg1(1,1,ind_k),
     &               nplw_gw,0.d0,tmp,npw)
          call dgemm('n','n',npw,npw,n2,1.d0,tmp,npw,s(n0_b+1,1),
     &               nplw_gw,0.d0,tmp1,npw)
          deallocate(tmp)
        endif
        if(nproc_b/=1) call dgop(tmp1,npw*npw,'  +',comm_b)
        ij=0
        do j=1,npw
          do i=1,j
            ij=ij+1
            if(ij>n1.and.ij<=n3) then
              p(ij-n1,ind_tau,ind_k)=p(ij-n1,ind_tau,ind_k)+tmp1(i,j)
            endif
          enddo
        enddo
        deallocate(tmp1)
      enddo   !! over ind_k
      deallocate(s)
      end


      subroutine p_ii_q_sclp_4r(pg1,nd_b,n0_b,ind_tau,pwsclp)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nd_b,n0_b,ind_tau
      real*8, intent(in) :: pg1(nplw_gw,nd_b,ndim3_k(me_k+1))
      real*8, intent(inout) :: pwsclp(nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: k,npw,i,ind_k,j,n2,n1,n3,ii,jj
      real*8, allocatable :: tmp(:,:),tmp1(:,:),s(:,:)
      n1=n0_b_pbint(me_b+1)
      n3=n1+nd_b_pbint(me_b+1)
      allocate(s(nplw_gw,nplw_gw))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        npw=nplwgw(k)
c ------ We form the overlap matrix --------------------------------
        call s_overlap_r(0,s,nplw_gw,npw,indpw_gw(1,k))
        n2=min(nd_b,npw-n0_b)
        allocate(tmp1(npw,npw))
        tmp1=0.d0
        if(n2>0) then
          allocate(tmp(npw,n2))
          call dgemm('t','n',npw,n2,npw,1.d0,s,nplw_gw,pg1(1,1,ind_k),
     &               nplw_gw,0.d0,tmp,npw)
          call dgemm('n','n',npw,npw,n2,1.d0,tmp,npw,s(n0_b+1,1),
     &               nplw_gw,0.d0,tmp1,npw)
          deallocate(tmp)
        endif
        if(nproc_b/=1) call dgop(tmp1,npw*npw,'  +',comm_b)
        do j=1,ncol_pw(k)
          jj=refcol_pw(j,k)-n_pbmt
          if(jj<1) cycle
          do i=1,nrow_pw(k)
            ii=refrow_pw(i,k)-n_pbmt
            if(ii<1) cycle
            pwsclp(i,j,ind_tau,ind_k)=pwsclp(i,j,ind_tau,ind_k)
     &                               +tmp1(ii,jj)
          enddo
        enddo
        deallocate(tmp1)
      enddo   !! over ind_k
      deallocate(s)
      end
