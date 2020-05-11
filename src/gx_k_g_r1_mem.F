      subroutine gx_k_g_r1_mem(g_c,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ispin
      complex*16, intent(out) :: g_c(nbasmpw,ndim3_k(me_k+1),nr_full,
     &                               nrel)
      integer :: k,nbnd,ind_k,nbas,irl,ibas0,i,i0,ind,j
      real*8 :: v(3),fc,pi2,de,gtm,green0
      complex*16 :: cc
      real*8, allocatable :: gx_tau(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp0(:,:),
     &                           tmp2(:,:)
      allocate(tmp0(nbndf,nbasmpw),tmp2(nbndf,nbndf))
      allocate(tmp1(nr_full,nbasmpw))
      allocate(tmp(nbasmpw,nbasmpw))
      pi2=pi+pi
      allocate(gx_tau(nbndf))
      g_c=(0.d0,0.d0)
c ---------------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        k=n3_mpi_k(me_k+1)+ind_k
        nbnd=n_bnd(k,ispin)
        if(ubi=='dft'.or.ubi==' hf') then
          do i=1,nbnd
            de=e_bnd(i,k,ispin)-chem_pot
            gtm=-green0(de,betta_t)
            do j=1,nbask(k)
              tmp0(i,j)=gtm*conjg(ev_bnd(j,i,ind_k,ispin))
            enddo
          enddo
        else
          call zgemm('n','c',nbnd,nbask(k),nbnd,(1.d0,0.d0),
     &               g_full_0(1,1,ind_k,ispin),nbndf,
     &               ev_bnd(1,1,ind_k,ispin),nbasmpw,(0.d0,0.d0),tmp0,
     &               nbndf)
        endif
        call zgemm('n','n',nbask(k),nbask(k),nbnd,(1.d0,0.d0),
     &             ev_bnd(1,1,ind_k,ispin),nbasmpw,tmp0,nbndf,
     &             (0.d0,0.d0),tmp,nbasmpw)
        nbas=nbask(k)/nrel
        do irl=1,nrel
          ibas0=(irl-1)*nbas
c ----- Fill in the first zone for second index -----------------------
          tmp1=(0.d0,0.d0)
          do i=1,nbas
            i0=indgb(i,k)
            call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
            tmp1(ind,1:nbask(k))=tmp(1:nbask(k),ibas0+i)
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),nbask(k),tmp1,-1)
c ----- Multiplication with the phase factor ------------------
          do i=1,nr_full
            v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &                               +r_pnt(3,i)*r0b0(:,3)
            fc=-pi2*dot_product(pnt(:,k),v)
            cc=exp(dcmplx(0.d0,fc))/amega
            g_c(1:nbask(k),ind_k,i,irl)=tmp1(i,1:nbask(k))*cc
          enddo
        enddo   !! over irl
      enddo   !! over ind_k
      deallocate(tmp0,tmp2,tmp,tmp1,gx_tau)
      end
