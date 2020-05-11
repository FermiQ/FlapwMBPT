      subroutine w_im_r_from_k(ind_tau)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      integer :: k,k0,ind,i0,j0,i,ind_k,ig,j,j1,iatom,isort,ia,ib,npw,
     &           ind_ir,ir,kst,ind_i,ii,jj
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: v_tmp(:,:,:),v_tmp1(:,:,:,:),tmp1(:,:),
     &                           tmp2(:,:),tmp3(:,:)
      allocate(v_tmp(nr_full,n_pbmt,ndim3_k(me_k+1)))
      allocate(tmp1(nr_full,n_pbmt))
      v_tmp=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        npw=nplwgw(k)
        if(npw/=0) then
          allocate(tmp2(n_pbmt,npw))
          if(ncmpl==1) then
            if(lsclp) then
              call p_c_from_r_mi_sclp(npw,tmp2,
     &                                pw_sclp(1,1,ind_tau,ind_k),k)
            else
              call p_c_from_r_mi(npw,tmp2,
     &                           p_w_q_tau_mi(1,1,ind_tau,ind_k),k)
            endif
          else
            tmp2=(0.d0,0.d0)
            if(lsclp) then
              do j=1,ncol_pw(k)
                jj=refcol_pw(j,k)-n_pbmt
                if(jj<=0) cycle
                do i=1,nrow_pw(k)
                  ii=refrow_pw(i,k)
                  if(ii>n_pbmt) cycle
                  call zcopy(1,pw_sclp(2*i-1,j,ind_tau,ind_k),1,
     &                       tmp2(ii,jj),1)
                enddo
              enddo
            else
              do i=1,npw
                call zcopy(nd_b_pbmt(me_b+1),
     &                     p_w_q_tau_mi(1,i,ind_tau,ind_k),1,
     &                     tmp2(n0_b_pbmt(me_b+1)+1,i),1)
              enddo
            endif
            if(nproc_b/=1) call dgop(tmp2,2*n_pbmt*npw,'  +',comm_b)
          endif
          tmp2=tmp2/2
          tmp1=(0.d0,0.d0)
          do i=1,npw
            i0=indpw_gw(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(ind,1:n_pbmt)=tmp2(1:n_pbmt,i)
          enddo
          deallocate(tmp2)
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp1,-1)
c ----- Multiplication with the phase factor ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
            fc=-2.d0*pi*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
            v_tmp(i,1:n_pbmt,ind_k)=tmp1(i,1:n_pbmt)*cc
          enddo
        endif
      enddo   !! over ind_k
      deallocate(tmp1)
c ----------------------------------------------------------------
      allocate(v_tmp1(k_star_max,ndim3_k(me_k+1),nr_red,n_pbmt))
      pi2=pi+pi
      v_tmp1=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me_k+1)
        k0=n3_mpi_k(me_k+1)+ind_k
        do kst=1,k_star(k0)
          k=k_list(kst,k0)
          ig=k_group(k)
          do j=1,nr_red
            j1=r_pnt_ind(j)
            j0=ip_r(j1,ig)
            call sym_vmi_left(ig,v_tmp(j0,:,ind_k),
     &                        v_tmp1(kst,ind_k,j,:))
            do iatom=1,natom
              isort=is(iatom)
              v=tshift(:,iatom,ig)-tshift_r(:,j1,ig)
              fc=pi2*dot_product(pnt(:,k),v)
              cc=dcmplx(cos(fc),sin(fc))
              ia=iopb(iatom)
              ib=ia+n_pbmt0(isort)-1
              v_tmp1(kst,ind_k,j,ia:ib)=v_tmp1(kst,ind_k,j,ia:ib)*cc
            enddo
          enddo
        enddo
      enddo
      deallocate(v_tmp)
c --------------------------------------------------------------------
      allocate(v_tmp(nqdiv,ndim3_nr_red(me_k+1),n_pbmt))
      allocate(tmp3(nqdiv,nr_red))
      do j=1,n_pbmt
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
      call fft3(ndiv(1),ndiv(2),ndiv(3),
     &          n_pbmt*ndim3_nr_red(me_k+1),v_tmp,1)
      v_tmp=v_tmp/dfloat(nqdiv)
c -----------------------------------------------------------------
      do j=1,n_pbmt
        tmp3=(0.d0,0.d0)
        do ind_i=1,ndim3_nr_red(me_k+1)
          i=n3_mpi_nr_red(me_k+1)+ind_i
          do ir=1,nqdiv
            tmp3(ir,i)=v_tmp(ir,ind_i,j)
          enddo
        enddo
        if(nproc_k/=1) call DGOP(tmp3,2*nr_red*nqdiv,'  +',comm_k)
        do ind_ir=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_ir
          w_rs_mi(:,j,ind_ir)=tmp3(ir,:)
        enddo
      enddo
      deallocate(v_tmp,tmp3)
      end
