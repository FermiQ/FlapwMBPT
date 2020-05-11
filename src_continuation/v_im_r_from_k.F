      subroutine v_im_r_from_k
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: k,k0,ind,i0,j0,i,ind_k,ig,j,j1,iatom,isort,ia,ib,
     &           ind_ir,ir,kst,ind_i,khbr,ind_j,npw,nd,n1,n2
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      real*8, allocatable :: tmp2r(:,:)
      complex*16, allocatable :: v_tmp(:,:,:),v_tmp1(:,:,:,:),tmp1(:,:),
     &	                         tmp2(:,:),tmp3(:,:),t1(:),t2(:),
     &                           sit(:,:,:)
      khbr=0
      if(ubi==' hf'.and.hybrid_type>2) khbr=1
      nd=ndim_tb_n_pbmt(me_t_b+1)
      n1=nmpi_tb_n_pbmt(me_t_b+1)+1
      n2=n1+nd-1
      if(nd/=0) then
        allocate(v_tmp(nr_full,nd,ndim3_k(me_k+1)))
        allocate(tmp1(nr_full,nd))
        v_tmp=(0.d0,0.d0)
      endif
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        npw=nplwgw(k)
        if(npw/=0) then
          if(nd/=0) allocate(tmp2(nd,npw))
          if(ncmpl==1) then
            allocate(tmp2r(n_pbmt,npw))
            tmp2r=0.d0
            if(nd/=0) then
              do j=1,npw
                if(khbr==0) call dcopy(nd,v_mi_3(1,j,ind_k),1,
     &                                 tmp2r(n1,j),1)
                if(khbr==1) call dcopy(nd,v_mi_3_sr(1,j,ind_k),1,
     &                                 tmp2r(n1,j),1)
              enddo
            endif
            if(nproc_t*nproc_b/=1) call dgop(tmp2r,n_pbmt*npw,'  +',
     &                                       comm_t_b)
            if(nd/=0) call v_c_from_r_mi_par(npw,tmp2,tmp2r,k,nd,n1,n2)
            deallocate(tmp2r)
          else
            if(nd/=0) then
              do j=1,npw
                if(khbr==0) call zcopy(nd,v_mi_3(1,j,ind_k),1,tmp2(1,j),
     &                                 1)
                if(khbr==1) call zcopy(nd,v_mi_3_sr(1,j,ind_k),1,
     &                                 tmp2(1,j),1)
              enddo
            endif
          endif
          if(nd/=0) then
            tmp1=(0.d0,0.d0)
            do i=1,npw
              i0=indpw_gw(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              tmp1(ind,:)=tmp2(:,i)
            enddo
            deallocate(tmp2)
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nd,tmp1,-1)
c ----- Multiplication with the phase factor ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		                           +r_pnt(3,i)*r0b0(:,3)
              fc=-2.d0*pi*dot_product(pnt(:,k),v)
              cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
              v_tmp(i,:,ind_k)=tmp1(i,:)*cc
            enddo
          endif
        endif
      enddo   !! over ind_k
      if(nd/=0) deallocate(tmp1)
c ------------------------------------------------------------------
      pi2=pi+pi
      if(nd/=0) then
        allocate(v_tmp1(k_star_max,ndim3_k(me_k+1),nr_red,nd))
        v_tmp1=(0.d0,0.d0)
      endif
      allocate(t1(n_pbmt))
      allocate(t2(n_pbmt))
      do ind_k=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_k
        do kst=1,k_star(k0)
          k=k_list(kst,k0)
          ig=k_group(k)
          do j=1,nr_red
            j1=r_pnt_ind(j)
            j0=ip_r(j1,ig)
            t1=(0.d0,0.d0)
            if(nd/=0) then
              do ind_i=1,nd
                i=nmpi_tb_n_pbmt(me_t_b+1)+ind_i
                t1(i)=v_tmp(j0,ind_i,ind_k)
              enddo
            endif
            if(nproc_t*nproc_b/=1) call dgop(t1,2*n_pbmt,'  +',comm_t_b)
            call sym_vmi_left(ig,t1,t2)
            do iatom=1,natom
              isort=is(iatom)
              v=tshift(:,iatom,ig)-tshift_r(:,j1,ig)
              fc=pi2*dot_product(pnt(:,k),v)
              cc=dcmplx(cos(fc),sin(fc))
              ia=iopb(iatom)
              ib=ia+n_pbmt0(isort)-1
              t2(ia:ib)=t2(ia:ib)*cc
            enddo
            if(nd/=0) then
              do ind_i=1,nd
                i=nmpi_tb_n_pbmt(me_t_b+1)+ind_i
                v_tmp1(kst,ind_k,j,ind_i)=t2(i)
              enddo
            endif
          enddo
        enddo
      enddo
      deallocate(t1,t2)
      if(nd/=0) deallocate(v_tmp)
c ------------------------------------------------------------------
      if(nd/=0) then
        allocate(v_tmp(nqdiv,ndim3_nr_red(me_k+1),nd))
        allocate(tmp3(nqdiv,nr_red))
        do j=1,nd
          tmp3=(0.d0,0.d0)
          do ind_k=1,ndim3_k(me_k+1)
            k0=n3_mpi_k(me_k+1)+ind_k
            do kst=1,k_star(k0)
              k=k_list(kst,k0)
              ig=k_group(k)
              call zone1_number(pnt(:,k),rb0,ndiv,ind)
              do i=1,nr_red
                tmp3(ind,i)=v_tmp1(kst,ind_k,i,j)
              enddo
            enddo
          enddo
          if(nproc_k/=1) call DGOP(tmp3,2*nr_red*nqdiv,'  +',comm_k)
          do ind_i=1,ndim3_nr_red(me_k+1)
            i=n3_mpi_nr_red(me_k+1)+ind_i
            do k=1,nqdiv
              v_tmp(k,ind_i,j)=tmp3(k,i)
            enddo
          enddo
        enddo
        deallocate(v_tmp1)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nd*ndim3_nr_red(me_k+1),v_tmp,
     &            1)
        v_tmp=v_tmp/dfloat(nqdiv)
c -----------------------------------------------------------------
        allocate(sit(nr_red,nd,ndim3_kk(me_k+1)))
        do ind_j=1,nd
          tmp3=(0.d0,0.d0)
          do ind_i=1,ndim3_nr_red(me_k+1)
            i=n3_mpi_nr_red(me_k+1)+ind_i
            do ir=1,nqdiv
              tmp3(ir,i)=v_tmp(ir,ind_i,ind_j)
            enddo
          enddo
          if(nproc_k/=1) call DGOP(tmp3,2*nr_red*nqdiv,'  +',comm_k)
          do ind_ir=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_ir
            sit(:,ind_j,ind_ir)=tmp3(ir,:)
          enddo
        enddo
        deallocate(v_tmp,tmp3)
      endif
c -----------------------------------------------------------------
      if(ndim_tb_nr_red(me_t_b+1)/=0) then
        allocate(s_it_r0(ndim_tb_nr_red(me_t_b+1),n_pbmt,
     &                   ndim3_kk(me_k+1)))
      endif
      allocate(tmp3(nr_red,n_pbmt))
      do ind_ir=1,ndim3_kk(me_k+1)
        tmp3=(0.d0,0.d0)
        if(nd/=0) then
          do ind_i=1,nd
            i=nmpi_tb_n_pbmt(me_t_b+1)+ind_i
            tmp3(:,i)=sit(:,ind_i,ind_ir)
          enddo
        endif
        if(nproc_t*nproc_b/=1) call DGOP(tmp3,2*nr_red*n_pbmt,'  +',
     &                                   comm_t_b)
        if(ndim_tb_nr_red(me_t_b+1)/=0) then
          do ind_i=1,ndim_tb_nr_red(me_t_b+1)
            i=nmpi_tb_nr_red(me_t_b+1)+ind_i
            s_it_r0(ind_i,:,ind_ir)=tmp3(i,:)
          enddo
        endif
      enddo
      deallocate(tmp3)
      if(nd/=0) deallocate(sit)
      end
