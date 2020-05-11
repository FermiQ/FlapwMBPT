      subroutine sigc_im_r_to_k(ispin,it,sig_mi,sig_tmp)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin,it
      complex*16, intent(in) :: sig_mi(nr_red,nfun,nrel,
     &	                               ndim3_kk(me_k+1))
      complex*16, intent(inout) :: sig_tmp(nbndf,nbndf,2,
     &                                     ndim3_k(me_k+1))
      integer :: ind,i,i0,j,ia,ib,ind_k,k,ig,jatom,jsort,nbnd,nbas,
     &	         irl,ind_ir,ir,k0,k1
      real*8 :: v(3),fc,pi2
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:),tmp4(:),
     &	                         p_c(:,:),pcc(:,:,:,:,:),pmi(:,:)
      allocate(p_c(nqdiv,nfun))
      allocate(pcc(nfun,nr_red,nrel,k_star_max,ndim3_k(me_k+1)))
      pcc=(0.d0,0.d0)
      do irl=1,nrel
        do i0=1,nr_red
          p_c=(0.d0,0.d0)
          do ind_ir=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_ir
            do j=1,nfun
              p_c(ir,j)=sig_mi(i0,j,irl,ind_ir)
            enddo
          enddo
          if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun,'  +',comm_k)
          call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,p_c,-1)
          do ind_k=1,ndim3_k(me_k+1)
            k=n3_mpi_k(me_k+1)+ind_k
            do k0=1,k_star(k)
              k1=k_list(k0,k)
              call zone1_number(pnt(1,k1),rb0,ndiv,ind)
              pcc(:,i0,irl,k0,ind_k)=pcc(:,i0,irl,k0,ind_k)+p_c(ind,:)
            enddo
          enddo
        enddo
      enddo
      deallocate(p_c)
      allocate(tmp(nr_full,nfun),tmp4(nfun))
      allocate(tmp1(nfun,nrel*nplw_gw))
      allocate(tmp2(nbndf,nbndf))
      allocate(pmi(nfun,nbndf))
      pi2=pi+pi
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbas=nbask(k)/nrel
        nbnd=n_bnd(k,ispin)
        pmi=(0.d0,0.d0)
        tmp1=(0.d0,0.d0)
        do irl=1,nrel
c ------- We fill in the first index completely -------------
          do i=1,nr_full
            i0=ir_ref(i)
            ig=r_group(i)
            ind=k_sym_0(k,ig)  ! IND=A^(-1)*K  in NQDIV-list
            ind=ind_k_star_back(ind) ! position in a given K-star
            if(irel<=1) then
              tmp4=pcc(:,i0,1,ind,ind_k)
            else if(irl==1) then
              tmp4=uj(1,ig)*pcc(:,i0,1,ind,ind_k)
     &            +uj(3,ig)*pcc(:,i0,2,ind,ind_k)
            else if(irl==2) then
              tmp4=uj(2,ig)*pcc(:,i0,1,ind,ind_k)
     &            +uj(4,ig)*pcc(:,i0,2,ind,ind_k)
            endif
            call sym_val_left(ig,tmp4,tmp(i,1:nfun),1)
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
c ----- Multiplication with the phase factors ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                               +r_pnt(3,i)*r0b0(:,3)
            fc=2.d0*pi*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))
            tmp(i,1:nfun)=tmp(i,1:nfun)*cc
          enddo
c -----------------------------------------------------------
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nfun,tmp,1)
          tmp=tmp/nr_full*sqrt(amega)
          do i=1,nplwgw(k)
            i0=indpw_gw(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            i0=(irl-1)*nplwgw(k)+i
            tmp1(1:nfun,i0)=tmp(ind,1:nfun)
          enddo
        enddo   !! over irl
        do irl=1,nrel
          i0=(irl-1)*nplwgw(k)+1
          call zgemm('n','n',nfun,nbnd,nplwgw(k),(1.d0,0.d0),
     &               tmp1(1,i0),nfun,sa(1,1,ind_k,irl,ispin),nplw_gw,
     &               (1.d0,0.d0),pmi,nfun)
        enddo   !! over irl
c -------- Contribution to Sigma ---------------------------
        call zgemm('c','n',nbnd,nbnd,nfun,(1.d0,0.d0),
     &	           z_bnd(1,1,ind_k,ispin),nfun,pmi,nfun,(0.d0,0.d0),
     &             tmp2,nbndf)
        do j=1,nbnd
          do i=1,j
            tmp2(i,j)=tmp2(i,j)+conjg(tmp2(j,i))
            if(i/=j) tmp2(j,i)=conjg(tmp2(i,j))
          enddo
        enddo
        sig_tmp(:,:,it,ind_k)=sig_tmp(:,:,it,ind_k)+tmp2
      enddo  !! over ind_k
      deallocate(tmp,tmp1,tmp2,tmp4,pcc,pmi)
      end
