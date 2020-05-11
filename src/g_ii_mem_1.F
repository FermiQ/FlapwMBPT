      subroutine g_ii_mem_1(it,ind_tau,g_c,ispin,nnn)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau,it,ispin,nnn
      complex*16, intent(out) :: g_c(nbasmpw,nnn,nrel,npnt)
      integer :: k,nbnd,ind_k,nbas,i,ind,i_tau,j,ib,ind_ista,ista,ii,
     &           ind_i,j0,jbas0,jrl,npr
      real*8 :: v(3),fc,pi2,tu
      complex*16 :: cc
      integer, allocatable :: nd_nb(:),n0_nb(:),nd_nb1(:),n0_nb1(:)
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp0(:,:),tmp2(:,:),gc(:,:,:),
     &                           gc1(:,:,:)
      allocate(nd_nb(nproc_b))
      allocate(n0_nb(nproc_b))
      call size_shift_par(nbasmpw,nproc_b,nd_nb,n0_nb)
      pi2=pi+pi
      allocate(gc(nd_nb(me_b+1),nbasmpw,ndim3_k(me_k+1)))
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbnd=n_bnd(k,ispin)
        allocate(tmp0(nd_nb(me_b+1),nbnd))
        tmp0=(0.d0,0.d0)
        if(ubi=='dft'.or.ubi==' hf') then
          i_tau=me_t*ndim3_tau+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          tu=tau_mesh(i_tau)
          allocate(gx_tau(nbnd))
          call g_x_tau(ispin,k,gx_tau,tu,nbnd,chem_pot)
          do ind_i=1,nd_nb(me_b+1)
            i=n0_nb(me_b+1)+ind_i
            if(i>nbask(k)) cycle
            do ib=1,nbnd
              tmp0(ind_i,ib)=ev_bnd(i,ib,ind_k,ispin)*gx_tau(ib)
            enddo
          enddo
          deallocate(gx_tau)
        else
          allocate(tmp2(nbnd,nbnd))
          call ferm_unpack_tau(tmp2,g_full(1,1,1,ind_tau,ind_k,ispin),
     &                         nbnd,nbnd,nbndf,it)
          call zgemm('n','n',nd_nb(me_b+1),nbnd,nbnd,(1.d0,0.d0),
     &               ev_bnd(n0_nb(me_b+1)+1,1,ind_k,ispin),nbasmpw,tmp2,
     &               nbnd,(0.d0,0.d0),tmp0,nd_nb(me_b+1))
          deallocate(tmp2)
        endif
        call zgemm('n','c',nd_nb(me_b+1),nbask(k),nbnd,(1.d0,0.d0),tmp0,
     &             nd_nb(me_b+1),ev_bnd(1,1,ind_k,ispin),nbasmpw,
     &             (0.d0,0.d0),gc(1,1,ind_k),nd_nb(me_b+1))
        deallocate(tmp0)
      enddo   !! over ind_k
c ------------------------------------------------------------------
      npr=nproc_b*nproc_k
      allocate(nd_nb1(npr))
      allocate(n0_nb1(npr))
      call size_shift_par(nbasmpw,npr,nd_nb1,n0_nb1)
      allocate(gc1(nd_nb1(me_k_b+1),nbasmpw,npnt))
      allocate(tmp0(nbasmpw,npnt))
      do j=1,nbasmpw
        tmp0=(0.d0,0.d0)
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          do ind_i=1,nd_nb(me_b+1)
            i=n0_nb(me_b+1)+ind_i
            tmp0(i,k)=gc(ind_i,j,ind_k)
          enddo
        enddo
        if(nproc_k/=1) call dgop(tmp0,2*nbasmpw*npnt,'  +',comm_k)
        if(nproc_b/=1) call dgop(tmp0,2*nbasmpw*npnt,'  +',comm_b)
        do k=1,npnt
          do ind_i=1,nd_nb1(me_k_b+1)
            i=n0_nb1(me_k_b+1)+ind_i
            gc1(ind_i,j,k)=tmp0(i,k)
          enddo
        enddo
      enddo
      deallocate(gc,tmp0,nd_nb,n0_nb)
c ------------------------------------------------------------------
      allocate(tmp0(nr_full,nd_nb1(me_k_b+1)))
      do k=1,npnt
        nbas=nbask(k)/nrel
        do jrl=1,nrel
          jbas0=(jrl-1)*nbas
c ----- Fill in the first zone for second index -----------------------
          tmp0=(0.d0,0.d0)
          do j=1,nbas
            j0=indgb(j,k)
            call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
            do ind_i=1,nd_nb1(me_k_b+1)
              tmp0(ind,ind_i)=gc1(ind_i,jbas0+j,k)
            enddo
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nd_nb1(me_k_b+1),tmp0,-1)
c ----- Multiplication with the phase factor -----------------------
          allocate(tmp2(nbasmpw,nr_full))
          tmp2=(0.d0,0.d0)
          do j=1,nr_full
            v(:)=r_pnt(1,j)*r0b0(:,1)+r_pnt(2,j)*r0b0(:,2)
     &                               +r_pnt(3,j)*r0b0(:,3)
            fc=-pi2*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))/sqrt(amega)
            do ind_i=1,nd_nb1(me_k_b+1)
              i=n0_nb1(me_k_b+1)+ind_i
              tmp2(i,j)=tmp0(j,ind_i)*cc
            enddo
          enddo
          if(npr/=1) call dgop(tmp2,2*nr_full*nbasmpw,'  +',comm_k_b)
          ii=0
          do ind_ista=1,ndim4_nr_red(me_k_b+1)
            ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
            do i=1,num_r_star(ista)
              ii=ii+1
              ind=ind_r_star(i,ista)
              g_c(:,ii,jrl,k)=tmp2(:,ind)
            enddo
          enddo
          deallocate(tmp2)
        enddo   !! over irl
      enddo   !! over ind_k
      deallocate(tmp0,gc1,nd_nb1,n0_nb1)
      end
