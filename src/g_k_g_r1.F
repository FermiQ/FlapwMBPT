      subroutine g_k_g_r1(it,ind_tau,g_c,ispin,nnn)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau,it,ispin,nnn
      complex*16, intent(out) :: g_c(nnn,npnt,max_r_inv,
     &                               ndim4_nr_red(me_k_b+1),nrel)
      integer :: k,nbnd,ind_k,nbas,irl,ibas0,i,i0,ind,i_tau,j,ib,
     &           ind_ista,ista,n0,n1,n2
      real*8 :: v(3),fc,pi2,tu
      complex*16 :: cc
      integer, allocatable :: nd_b(:),n0_b(:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp0(:,:),
     &                           tmp2(:,:),gc(:,:,:,:)
      pi2=pi+pi
      allocate(nd_b(nproc_b))
      allocate(n0_b(nproc_b))
      call size_shift_par(nbasmpw,nproc_b,nd_b,n0_b)
      n0=n0_b(me_b+1)
      n1=n0+nd_b(me_b+1)
      allocate(gc(nd_b(me_b+1),ndim3_k(me_k+1),nr_full_inv,nrel))
      gc=(0.d0,0.d0)
c ---------------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbnd=n_bnd(k,ispin)
        if(n0>=nbask(k)) cycle
        n2=min(nd_b(me_b+1),nbask(k)-n0)
        allocate(tmp0(nbnd,nbask(k)))
        if(ubi=='dft'.or.ubi==' hf') then
          i_tau=me_t*ndim3_tau+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          tu=tau_mesh(i_tau)
          allocate(gx_tau(nbnd))
          call g_x_tau(ispin,k,gx_tau,tu,nbnd,chem_pot)
          do j=1,nbask(k)
            do ib=1,nbnd
              tmp0(ib,j)=gx_tau(ib)*conjg(ev_bnd(j,ib,ind_k,ispin))
            enddo
          enddo
          deallocate(gx_tau)
        else
          allocate(tmp2(nbnd,nbnd))
          call ferm_unpack_tau(tmp2,g_full(1,1,1,ind_tau,ind_k,ispin),
     &                         nbnd,nbnd,nbndf,it)
          call zgemm('n','c',nbnd,nbask(k),nbnd,(1.d0,0.d0),tmp2,nbnd,
     &               ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),tmp0,
     &               nbnd)
          deallocate(tmp2)
        endif
        allocate(tmp(n2,nbasmpw))
        call zgemm('n','n',n2,nbask(k),nbnd,(1.d0,0.d0),
     &             ev_bnd(n0+1,1,ind_k,ispin),nbasmpw,tmp0,nbnd,
     &             (0.d0,0.d0),tmp,n2)
        deallocate(tmp0)
        nbas=nbask(k)/nrel
        allocate(tmp1(nr_full,n2))
        do irl=1,nrel
          ibas0=(irl-1)*nbas
c ----- Fill in the first zone for second index -----------------------
          tmp1=(0.d0,0.d0)
          do i=1,nbas
            i0=indgb(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(ind,1:n2)=tmp(1:n2,ibas0+i)
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),n2,tmp1,-1)
c ----- Multiplication with the phase factor ------------------
          do i0=1,nr_full_inv
            i=i_r_inv(i0)
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                               +r_pnt(3,i)*r0b0(:,3)
            fc=-pi2*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))/amega
            gc(1:n2,ind_k,i0,irl)=tmp1(i,1:n2)*cc
          enddo
        enddo   !! over irl
        deallocate(tmp,tmp1)
      enddo   !! over ind_k
      allocate(tmp(npnt,nr_full_inv))
      do irl=1,nrel
        do i=1,nbasmpw
          tmp=(0.d0,0.d0)
          if(i>n0.and.i<=n1) then
            do ind_k=1,ndim3_k(me_k+1)
              k=n3_mpi_k(me_k+1)+ind_k
              tmp(k,:)=gc(i-n0,ind_k,:,irl)
            enddo
          endif
          if(nproc_k*nproc_b/=1) call dgop(tmp,2*nr_full_inv*npnt,'  +',
     &                                     comm_k_b)
          do ind_ista=1,ndim4_nr_red(me_k_b+1)
            ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
            do j=1,num_r_inv(ista)
              ind=ind_r_inv(j,ista)
              g_c(i,:,j,ind_ista,irl)=tmp(:,ind)
            enddo
          enddo
        enddo
      enddo
      deallocate(gc,tmp,n0_b,nd_b)
      end
