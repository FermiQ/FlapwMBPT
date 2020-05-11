      subroutine sigx_im_k_from_r(sig_mi)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      complex*16, intent(inout) :: sig_mi(ndim_tb_nr_red(me_t_b+1),nfun,
     &                                    nrel,ndim3_kk(me_k+1))
      integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,irl,nbas,
     &	         ibas0,ind_ir,ir,k0,k1,ic,j0,ind_i
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),s(:,:),tmp3(:),
     &	                         p_c(:,:),pcc(:,:,:,:,:),tmp2(:,:)
      if(ndim_tb_nr_red(me_t_b+1)/=0) then
        allocate(p_c(nqdiv,nfun))
        allocate(pcc(nfun,ndim_tb_nr_red(me_t_b+1),nrel,k_star_max,
     &               ndim3_k(me_k+1)))
        do irl=1,nrel
          do ind_i=1,ndim_tb_nr_red(me_t_b+1)
            p_c=(0.d0,0.d0)
            do ind_ir=1,ndim3_kk(me_k+1)
              ir=n3_mpi_kk(me_k+1)+ind_ir
              do j=1,nfun
                p_c(ir,j)=sig_mi(ind_i,j,irl,ind_ir)
              enddo
            enddo
            if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun,'  +',comm_k)
            call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,p_c,-1)
            do ind_k=1,ndim3_k(me_k+1)
              k=n3_mpi_k(me_k+1)+ind_k
              do k0=1,k_star(k)
                k1=k_list(k0,k)
                call zone1_number(pnt(1,k1),rb0,ndiv,ind)
                pcc(:,ind_i,irl,k0,ind_k)=p_c(ind,:)
              enddo
            enddo
          enddo
        enddo
        deallocate(p_c)
      endif
      pi2=pi+pi
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbas=nbask(k)/nrel
        do irl=1,nrel
          ibas0=nbas*(irl-1)
c ------- We fill in the first index completely -------------
          allocate(tmp(nr_full,nfun),tmp3(nfun))
          tmp=(0.d0,0.d0)
          if(ndim_tb_nr_red(me_t_b+1)/=0) then
            do i=1,nr_full
              i0=ir_ref(i)
              ind_i=i0-nmpi_tb_nr_red(me_t_b+1)
              if(ind_i<=0.or.ind_i>ndim_tb_nr_red(me_t_b+1)) cycle
              ig=r_group(i)
              ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
              ind=ind_k_star_back(ind) ! position in a given K-star
              if(irel<=1) then
                tmp3=pcc(:,ind_i,1,ind,ind_k)
              else if(irl==1) then
                tmp3=uj(1,ig)*pcc(:,ind_i,1,ind,ind_k)
     &              +uj(3,ig)*pcc(:,ind_i,2,ind,ind_k)
              else if(irl==2) then
                tmp3=uj(2,ig)*pcc(:,ind_i,1,ind,ind_k)
     &              +uj(4,ig)*pcc(:,ind_i,2,ind,ind_k)
              endif
              call sym_val_left(ig,tmp3,tmp(i,:),1)
              do jatom=1,natom
                jsort=is(jatom)
                v=tshift(:,jatom,ig)-tshift_r(:,i,ig)
                fc=pi2*dot_product(pnt(:,k),v)
                cc=dcmplx(cos(fc),sin(fc))
                ia=io_lem(jatom)
                ib=ia+lfunm(jsort)-1
                tmp(i,ia:ib)=tmp(i,ia:ib)*cc
              enddo
            enddo
          endif
          deallocate(tmp3)
          if(nproc_t*nproc_b/=1) call dgop(tmp,2*nr_full*nfun,'  +',
     &                                     comm_t_b)
c ----- Multiplication with the phase factors ------------------
          if(ndim_tb_nfun(me_t_b+1)/=0) then
            allocate(tmp2(nr_full,ndim_tb_nfun(me_t_b+1)))
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &		                           +r_pnt(3,i)*r0b0(:,3)
              fc=pi2*dot_product(pnt(:,k),v)
              cc=exp(dcmplx(0.d0,fc))
              do ind_i=1,ndim_tb_nfun(me_t_b+1)
                i0=nmpi_tb_nfun(me_t_b+1)+ind_i
                tmp2(i,ind_i)=tmp(i,i0)*cc
              enddo
            enddo
          endif
          deallocate(tmp)
c -----------------------------------------------------------
          if(ndim_tb_nfun(me_t_b+1)/=0) then
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),ndim_tb_nfun(me_t_b+1),
     &                tmp2,1)
            tmp2=tmp2/nr_full*sqrt(amega)
            allocate(tmp1(ndim_tb_nfun(me_t_b+1),nplw_gw))
            tmp1=(0.d0,0.d0)
            do i=1,nplwgw(k)
              i0=indpw_gw(i,k)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              tmp1(:,i)=tmp2(ind,:)
            enddo
            deallocate(tmp2)
c ------ We form the overlap matrix -----------------------------
            allocate(s(nplw_gw,nbasmpw/nrel))
            do j=1,nbas
              j0=indgb(j,k)
              do i=1,nplwgw(k)
                i0=indpw_gw(i,k)
                ia=igbs(1,j0)-igbs(1,i0)
                ib=igbs(2,j0)-igbs(2,i0)
                ic=igbs(3,j0)-igbs(3,i0)
                ind=indplw(ia,ib,ic)
                if(complex_ro) then
                  s(i,j)=dcmplx(sovr(ind),sovi(ind))
                else
                  s(i,j)=dcmplx(sovr(ind),0.d0)
                endif
              enddo
            enddo
            call zgemm('n','n',ndim_tb_nfun(me_t_b+1),nbas,nplwgw(k),
     &                 (1.d0,0.d0),tmp1,ndim_tb_nfun(me_t_b+1),
     &                 s,nplw_gw,(0.d0,0.d0),pvg(1,ibas0+1,ind_k),
     &                 ndim_tb_nfun(me_t_b+1))
            deallocate(s,tmp1)
          endif
        enddo   !! over irl
      enddo   !! over ind_k
      if(ndim_tb_nr_red(me_t_b+1)/=0) deallocate(pcc)
      end
