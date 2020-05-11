      subroutine v_mm_from_q_sort(v_rs_mm,n1pbm,n2pbm,nd_nrr,n0_nrr,
     &                            n1pbmt,n2pbmt,nat1,nat2,isort,jsort,
     &                            nrr_reds)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n1pbm,n2pbm,nd_nrr,n1pbmt,n2pbmt,nat1,nat2,
     &                       isort,jsort,n0_nrr,
     &                       nrr_reds(3,nqdiv*nat1*nat2)
      real*8, intent(out) :: v_rs_mm(n1pbm,n2pbm,nd_nrr)
      integer :: k,k0,ind,iatom,i0,j0,ir0,ir,jatom,ig,ind_ir0,ind_i,j,i,
     &           ind_j,m0,m,jj0,jj,ind_q,len,khbr,ind0,jnd0,n0,iat,npr
      real*8 :: pi2,tt(3),phase
      integer, allocatable :: nd_pb(:),n0_pb(:)
      complex*16, allocatable :: w0(:,:),v1(:,:),tr1(:),tr2(:),
     &                           w1(:,:,:),w2(:,:)
      pi2=pi+pi
      khbr=0
      if(ubi==' hf'.and.hybrid_type>2) khbr=1
      v_rs_mm=0.d0
      npr=nproc_t*nproc_b
      allocate(nd_pb(npr))
      allocate(n0_pb(npr))
      call size_shift_par(n2pbmt,npr,nd_pb,n0_pb)
      if(nd_pb(me_t_b+1)/=0) allocate(w1(nqdiv,n1pbmt,nd_pb(me_t_b+1)))
      allocate(w0(n1pbmt,n2pbmt))
      allocate(w2(n1pbmt,n2pbmt))
      allocate(tr1(nat1))
      allocate(tr2(nat2))
      ind0=iopb(iat_1(isort))-1
      jnd0=iopb(iat_1(jsort))-1
      if(nd_pb(me_t_b+1)/=0) w1=(0.d0,0.d0)
      do ind_q=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_q
        if(khbr==0) then
          if(ncmpl==1) then
            call v_coul_mm_sort_r(w0,ind_q,v_mm_3(1,ind_q),n1pbmt,
     &                            n2pbmt,ind0,jnd0,isort,jsort,nat1,
     &                            nat2,n1pbm,n2pbm)
          else
            call v_coul_mm_sort(w0,ind_q,v_mm_3(1,ind_q),n1pbmt,n2pbmt,
     &                          ind0,jnd0)
          endif
        else if(khbr==1) then
          if(ncmpl==1) then
            call v_coul_mm_sort_r(w0,ind_q,v_mm_3_sr(1,ind_q),n1pbmt,
     &                            n2pbmt,ind0,jnd0,isort,jsort,nat1,
     &                            nat2,n1pbm,n2pbm)
          else
            call v_coul_mm_sort(w0,ind_q,v_mm_3_sr(1,ind_q),n1pbmt,
     &                          n2pbmt,ind0,jnd0)
          endif
        endif
        do m=1,k_star(k0)
          k=k_list(m,k0)
          ig=k_group(k)
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
          do iatom=1,nat1
            iat=iat_sort(iatom,isort)
            tt=tshift(:,iat,ig)
            phase=pi2*dot_product(pnt(:,k),tt)
            tr1(iatom)=dcmplx(cos(phase),sin(phase))
          enddo
          allocate(v1(n1pbmt,n2pbmt))
          v1=(0.d0,0.d0)
          if(nd_pb(me_t_b+1)/=0) then
            do ind_i=1,nd_pb(me_t_b+1)
              i=n0_pb(me_t_b+1)+ind_i
              v1(:,i)=w0(:,i)
            enddo
            call sym_w_mt_sort(1,ig,tr1,v1(1,n0_pb(me_t_b+1)+1),
     &                         nd_pb(me_t_b+1),nat1,n1pbmt,isort)
          endif
          if(npr/=1) call dgop(v1,2*n1pbmt*n2pbmt,'  +',comm_t_b)
          w2=v1
          deallocate(v1)
          allocate(v1(n2pbmt,n1pbmt))
          v1=transpose(w2)
          do iatom=1,nat2
            iat=iat_sort(iatom,jsort)
            tt=tshift(:,iat,ig)
            phase=pi2*dot_product(pnt(:,k),tt)
            tr2(iatom)=dcmplx(cos(phase),sin(phase))
          enddo
          call sym_w_mt_sort(2,ig,tr2,v1,n1pbmt,nat2,n2pbmt,jsort)
          w2=transpose(v1)
          deallocate(v1)
          if(nd_pb(me_t_b+1)/=0) then
            do ind_j=1,nd_pb(me_t_b+1)
              j=n0_pb(me_t_b+1)+ind_j
              w1(ind,:,ind_j)=w2(:,j)
            enddo
          endif
        enddo   !! over m
      enddo  !! over ind_q
      deallocate(tr1,tr2,w0,w2)
      if(nd_pb(me_t_b+1)/=0) then
        if(nproc_k/=1) call dgop(w1,2*nqdiv*n1pbmt*nd_pb(me_t_b+1),
     &                           '  +',comm_k)
c -----------------------------------------------------------------
        do ind_j=1,nd_pb(me_t_b+1)
          call fft3(ndiv(1),ndiv(2),ndiv(3),n1pbmt,w1(1,1,ind_j),1)
        enddo
        w1=w1/dfloat(nqdiv)
        m0=n0_pb(me_t_b+1)+1
        m=m0+nd_pb(me_t_b+1)-1
        do ind_ir0=1,nd_nrr
          ir0=n0_nrr+ind_ir0
          ir=nrr_reds(1,ir0)
          iatom=nrr_reds(2,ir0)
          jatom=nrr_reds(3,ir0)
          n0=n2pbm*(jatom-1)+1
          if(n0>m) cycle
          j=n0+n2pbm-1
          if(j<m0) cycle
          j0=max(n0,m0)
          j=min(j,m)
          jj0=j0-n0_pb(me_t_b+1)
          jj=j-n0_pb(me_t_b+1)
          j0=j0-n0+1
          j=j-n0+1
          i0=n1pbm*(iatom-1)+1
          i=i0+n1pbm-1
          v_rs_mm(:,j0:j,ind_ir0)=w1(ir,i0:i,jj0:jj)
        enddo
      endif
      len=n1pbm*n2pbm*nd_nrr
      if(npr/=1.and.len/=0) call dgop(v_rs_mm,len,'  +',comm_t_b)
      if(nd_pb(me_t_b+1)/=0) deallocate(w1)
      deallocate(nd_pb,n0_pb)
      end
