      subroutine gx_rs_mi_mem(gr,gc,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin
      real*8, intent(out) :: gr(ndim_tb_nr_red(me_t_b+1),nfun,
     &                          ndim3_kk(me_k+1))
      complex*16, intent(out) :: gc(ndim_tb_nr_red(me_t_b+1),nfun,
     &                              ndim3_kk(me_k+1),nrel)
      integer :: k,k0,nbnd,ind,i,i0,ind_k,ig,j,j1,j0,iatom,isort,ia,ib,
     &           i_star,nbas,irl,ii,ind_ir,ir,iat,ie,iii,ind_j,
     &           ns,indls,jjj,l,li,ls,m,nd,nee,nls
      real*8 :: v(3),fc,pi2,de,green0,gtm
      complex*16 :: cc
      integer, allocatable :: adr_ls(:,:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: g_tmp(:,:,:),tmp(:,:),tmp1(:,:),
     &                           tmp0(:,:),tmp2(:,:),tmp3(:),
     &                           g_mi(:,:,:,:),z_ord(:,:)
      pi2=pi+pi
      allocate(adr_ls(2,nsort*(nrel*maxb+1)))
      ind=0
      do isort=1,nsort
        do li=1,nrel*lmb(isort)+1
          ind=ind+1
          adr_ls(1,ind)=li
          adr_ls(2,ind)=isort
        enddo
      enddo
      nls=ind
      if(irel/=2) gr=0.d0
      if(irel==2) gc=(0.d0,0.d0)
c ------------------------------------------------------------------
      do ls=1,nls
        isort=adr_ls(2,ls)
        li=adr_ls(1,ls)
        if(irel<2) then
          l=li-1
          nd=l+l+1
        else if(irel==2) then
          l=li/2
          iii=2*(li-l-l+2)-5
          jjj=l+l+iii
          nd=((li+1)/2)*2
        endif
        nee=ntle(l,isort)
        ie=n_e_apw(l,isort)
        if(augm(ie,l,isort)=='LAW') nee=nee+1
        if(augm(ie,l,isort)=='AWL') nee=nee+1
        ns=nd*nee*nhsort(isort)
        indls=nfun_ls_adr(li,isort)
        if(ndim_tb_nr_red(me_t_b+1)/=0) then
          allocate(g_mi(nqdiv,ndim_tb_nr_red(me_t_b+1),ns,nrel))
        endif
        allocate(tmp(nr_full,ns),tmp1(ns,nbasmpw))
        allocate(g_tmp(nr_full,ns,nrel),tmp3(ns))
        allocate(tmp0(ns,nbndf),tmp2(nbndf,nbndf))
        if(ndim_tb_nr_red(me_t_b+1)/=0) g_mi=(0.d0,0.d0)
        allocate(gx_tau(nbndf))
        do ind_k=1,ndim3_k(me_k+1)
          k0=n3_mpi_k(me_k+1)+ind_k
          nbnd=n_bnd(k0,ispin)
          nbas=nbask(k0)/nrel
          allocate(z_ord(ns,nbndf))
          do j=1,ns
            j0=nfun_order(indls+j-1)
            do ib=1,nbnd
              z_ord(j,ib)=z_bnd(j0,ib,ind_k,ispin)
            enddo
          enddo
c -----------------------------------------------------------------
          g_tmp=(0.d0,0.d0)
          if(ubi=='dft'.or.ubi==' hf') then
            do j=1,nbnd
              de=e_bnd(j,k0,ispin)-chem_pot
              gtm=-green0(de,betta_t)
              do i=1,ns
                tmp0(i,j)=z_ord(i,j)*gtm
              enddo
            enddo
          else
            call zgemm('n','n',ns,nbnd,nbnd,(1.d0,0.d0),z_ord,ns,
     &                 g_full_0(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),tmp0,
     &                 ns)
          endif
          deallocate(z_ord)
          call zgemm('n','c',ns,nbask(k0),nbnd,(1.d0,0.d0),tmp0,ns,
     &               ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),tmp1,
     &               ns)
          tmp1=tmp1/sqrt(amega)
          do irl=1,nrel
            tmp=(0.d0,0.d0)
            do i=1,nbas
              i0=indgb(i,k0)
              call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
              ii=(irl-1)*nbas+i
              tmp(ind,1:ns)=tmp1(1:ns,ii)
            enddo
            call fft3(nrdiv(1),nrdiv(2),nrdiv(3),ns,tmp,-1)
c ----- Multiplication with the phase factor ------------------
            do i=1,nr_full
              v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                                 +r_pnt(3,i)*r0b0(:,3)
              fc=-pi2*dot_product(pnt(:,k0),v)
              cc=exp(dcmplx(0.d0,fc))
              g_tmp(i,1:ns,irl)=tmp(i,1:ns)*cc
            enddo
          enddo   !! over irl
c ----------------------------------------------------------------
          do i_star=1,k_star(k0)
            k=k_list(i_star,k0)
            ig=k_group(k)
            call zone1_number(pnt(:,k),rb0,ndiv,ind)
            if(ndim_tb_nr_red(me_t_b+1)/=0) then
              do ind_j=1,ndim_tb_nr_red(me_t_b+1)
                j=nmpi_tb_nr_red(me_t_b+1)+ind_j
                j1=r_pnt_ind(j)
                j0=ip_r(j1,ig)
                do irl=1,nrel
                  if(irel<=1) then
                    tmp3=g_tmp(j0,:,1)
                  else if(irl==1) then
                    tmp3=uj(1,ig)*g_tmp(j0,:,1)+uj(3,ig)*g_tmp(j0,:,2)
                  else if(irl==2) then
                    tmp3=uj(2,ig)*g_tmp(j0,:,1)+uj(4,ig)*g_tmp(j0,:,2)
                  endif
                  call sym_left_tl(ig,tmp3,g_mi(ind,ind_j,:,irl),1,li,
     &                             isort,ns)
                  do iat=1,nhsort(isort)
                    iatom=iat_sort(iat,isort)
                    v=tshift(:,iatom,ig)-tshift_r(:,j1,ig)
                    fc=pi2*dot_product(pnt(:,k),v)
                    cc=dcmplx(cos(fc),sin(fc))
                    do ie=1,nee
                      do m=1,nd
                        ia=indback(m,ie,iat,li,isort)
                        g_mi(ind,ind_j,ia,irl)=g_mi(ind,ind_j,ia,irl)*cc
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo   !! over ind_k
        deallocate(g_tmp,tmp0,tmp,tmp1,tmp2,tmp3,gx_tau)
        if(ndim_tb_nr_red(me_t_b+1)/=0) then
          if(nproc_k/=1) call dgop(g_mi,
     &                         2*nqdiv*ndim_tb_nr_red(me_t_b+1)*ns*nrel,
     &                             '  +',comm_k)
          call fft3(ndiv(1),ndiv(2),ndiv(3),
     &              ns*ndim_tb_nr_red(me_t_b+1)*nrel,g_mi,1)
          g_mi=g_mi/dfloat(nqdiv)
          do ind_ir=1,ndim3_kk(me_k+1)
            ir=n3_mpi_kk(me_k+1)+ind_ir
            do j=1,ns
              j0=nfun_order(indls+j-1)
              if(irel<=1) then
                gr(:,j0,ind_ir)=g_mi(ir,:,j,1)
              else if(irel==2) then
                do irl=1,2
                  gc(:,j0,ind_ir,irl)=g_mi(ir,:,j,irl)
                enddo
              endif
            enddo
          enddo
          deallocate(g_mi)
        endif
      enddo   !! over ind_ls
      deallocate(adr_ls)
      end
