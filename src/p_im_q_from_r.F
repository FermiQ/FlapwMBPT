      subroutine p_im_q_from_r(p_mi,ind_tau,pmi)
c     key = 0 : contribution to P_W_SOLID_Q_TAU -----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      real*8, intent(in) :: p_mi(nr_red,n_pbmt,ndim3_kk(me_k+1))
      real*8, intent(inout) :: pmi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                             ndim3_tn,ndim3_k(me_k+1))
      integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,ind_ir,ir,k0,
     &           k1,npw,nd,n0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      real*8, allocatable :: tmp1r(:,:),sr(:,:)
      complex*16, allocatable :: tmp(:,:),p_c(:,:),pcc(:,:,:,:),
     &                           tmp1(:,:),s(:,:)
      nd=nd_b_pbmt(me_b+1)
      n0=n0_b_pbmt(me_b+1)
      pi2=pi+pi
      allocate(p_c(nqdiv,n_pbmt))
      allocate(pcc(n_pbmt,nr_red,k_star_max,ndim3_k(me_k+1)))
      do i0=1,nr_red
        p_c=(0.d0,0.d0)
        do ind_ir=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_ir
          do j=1,n_pbmt
            p_c(ir,j)=p_mi(i0,j,ind_ir)
          enddo
        enddo
        if(nproc_k/=1) call dgop(p_c,2*nqdiv*n_pbmt,'  +',comm_k)
        call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,p_c,-1)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do k0=1,k_star(k)
            k1=k_list(k0,k)
            call zone1_number(pnt(1,k1),rb0,ndiv,ind)
            pcc(:,i0,k0,ind_k)=p_c(ind,:)
          enddo
        enddo
      enddo  !! over i0
      deallocate(p_c)
c ---------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
c ------- We fill in the first index completely -------------
        allocate(tmp(nr_full,n_pbmt))
        do i=1,nr_full
          i0=ir_ref(i)
          ig=r_group(i)
          ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
          ind=ind_k_star_back(ind) ! position in a given K-star
          call sym_vmi_left(ig,pcc(:,i0,ind,ind_k),tmp(i,:))
          do jatom=1,natom
            jsort=is(jatom)
            v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
            fc=pi2*dot_product(pnt(:,k),v)
            cc=dcmplx(cos(fc),sin(fc))
            ia=iopb(jatom)
            ib=ia+n_pbmt0(jsort)-1
            tmp(i,ia:ib)=tmp(i,ia:ib)*cc
          enddo
c ----- Multiplication with the phase factors ------------------
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          fc=2.d0*pi*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))
          tmp(i,1:n_pbmt)=tmp(i,1:n_pbmt)*cc
        enddo   !! over ist
c -----------------------------------------------------------
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp,1)
        tmp=tmp/nr_full*sqrt(amega)
        npw=nplwgw(k)
        if(ncmpl==1) then
          allocate(tmp1r(n_pbmt,npw))
          call p_r_from_c_mi(npw,tmp,tmp1r,k)
          deallocate(tmp)
          allocate(sr(npw,npw))
          call s_overlap_r(0,sr,npw,npw,indpw_gw(1,k))
          call dgemm('n','n',nd,npw,npw,1.d0,tmp1r(n0+1,1),n_pbmt,sr,
     &               npw,1.d0,pmi(1,1,ind_tau,ind_k),nd)
          deallocate(sr,tmp1r)
        else
          allocate(tmp1(n_pbmt,npw))
          do i=1,npw
            i0=indpw_gw(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(1:n_pbmt,i)=tmp(ind,1:n_pbmt)
          enddo
          deallocate(tmp)
          allocate(s(npw,npw))
          call s_overlap(0,s,npw,npw,indpw_gw(1,k))
          s=s/amega
          call zgemm('n','n',nd,npw,npw,(1.d0,0.d0),tmp1(n0+1,1),n_pbmt,
     &               s,npw,(1.d0,0.d0),pmi(1,1,ind_tau,ind_k),nd)
          deallocate(s,tmp1)
        endif
      enddo  !! over ind_k
      deallocate(pcc)
      end


      subroutine p_im_q_sclp(p_mi,ind_tau,pwsclp)
c     key = 0 : contribution to P_W_SOLID_Q_TAU -----------------
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau
      real*8, intent(in) :: p_mi(nr_red,n_pbmt,ndim3_kk(me_k+1))
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,ind_ir,ir,k0,
     &           k1,npw,nd,n0
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      real*8, allocatable :: tmp1r(:,:),sr(:,:),s1(:,:)
      complex*16, allocatable :: tmp(:,:),p_c(:,:),pcc(:,:,:,:),
     &                           tmp1(:,:),s(:,:)
      nd=nd_b_pbmt(me_b+1)
      n0=n0_b_pbmt(me_b+1)
      pi2=pi+pi
      allocate(p_c(nqdiv,n_pbmt))
      allocate(pcc(n_pbmt,nr_red,k_star_max,ndim3_k(me_k+1)))
      do i0=1,nr_red
        p_c=(0.d0,0.d0)
        do ind_ir=1,ndim3_kk(me_k+1)
          ir=n3_mpi_kk(me_k+1)+ind_ir
          do j=1,n_pbmt
            p_c(ir,j)=p_mi(i0,j,ind_ir)
          enddo
        enddo
        if(nproc_k/=1) call dgop(p_c,2*nqdiv*n_pbmt,'  +',comm_k)
        call fft3(ndiv(1),ndiv(2),ndiv(3),n_pbmt,p_c,-1)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do k0=1,k_star(k)
            k1=k_list(k0,k)
            call zone1_number(pnt(1,k1),rb0,ndiv,ind)
            pcc(:,i0,k0,ind_k)=p_c(ind,:)
          enddo
        enddo
      enddo  !! over i0
      deallocate(p_c)
c ---------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
c ------- We fill in the first index completely -------------
        allocate(tmp(nr_full,n_pbmt))
        do i=1,nr_full
          i0=ir_ref(i)
          ig=r_group(i)
          ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
          ind=ind_k_star_back(ind) ! position in a given K-star
          call sym_vmi_left(ig,pcc(:,i0,ind,ind_k),tmp(i,:))
          do jatom=1,natom
            jsort=is(jatom)
            v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
            fc=pi2*dot_product(pnt(:,k),v)
            cc=dcmplx(cos(fc),sin(fc))
            ia=iopb(jatom)
            ib=ia+n_pbmt0(jsort)-1
            tmp(i,ia:ib)=tmp(i,ia:ib)*cc
          enddo
c ----- Multiplication with the phase factors ------------------
          v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                             +r_pnt(3,i)*r0b0(:,3)
          fc=2.d0*pi*dot_product(pnt(:,k),v)
          cc=exp(dcmplx(0.d0,fc))
          tmp(i,1:n_pbmt)=tmp(i,1:n_pbmt)*cc
        enddo   !! over i
c -----------------------------------------------------------
        call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n_pbmt,tmp,1)
        tmp=tmp/nr_full*sqrt(amega)
        npw=nplwgw(k)
        if(ncmpl==1) then
          allocate(tmp1r(n_pbmt,npw))
          call p_r_from_c_mi(npw,tmp,tmp1r,k)
          deallocate(tmp)
          allocate(sr(npw,npw))
          call s_overlap_r(0,sr,npw,npw,indpw_gw(1,k))
          allocate(s1(n_pbmt,npw))
          call dgemm('n','n',n_pbmt,npw,npw,1.d0,tmp1r,n_pbmt,sr,
     &               npw,0.d0,s1,n_pbmt)
          deallocate(sr,tmp1r)
          call p_r_add_mi_sclp(pwsclp(1,1,ind_tau,ind_k),s1,npw,k)
          deallocate(s1)
        else
          allocate(tmp1(n_pbmt,npw))
          do i=1,npw
            i0=indpw_gw(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(1:n_pbmt,i)=tmp(ind,1:n_pbmt)
          enddo
          deallocate(tmp)
          allocate(s(npw,npw))
          call s_overlap(0,s,npw,npw,indpw_gw(1,k))
          s=s/amega
          allocate(tmp(n_pbmt,npw))
          call zgemm('n','n',n_pbmt,npw,npw,(1.d0,0.d0),tmp1,n_pbmt,
     &               s,npw,(0.d0,0.d0),tmp,n_pbmt)
          deallocate(s,tmp1)
          call p_c_add_mi_sclp(pwsclp(1,1,ind_tau,ind_k),tmp,npw,k)
          deallocate(tmp)
        endif
      enddo  !! over ind_k
      deallocate(pcc)
      end
