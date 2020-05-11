      subroutine p_ii_q_mem_3r(pg,pg1,nd_b,n0_b,nnn)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: nnn,nd_b,n0_b
      complex*16, intent(in) :: pg(nplw_gw,nnn,npnt)
      real*8, intent(out) :: pg1(nplw_gw,nd_b,ndim3_k(me_k+1))
      integer :: k,npw,i,ii,ind,ind_ista,ind_j,ind_k,ista,j,j0,nsta
      real*8 :: v(3),pi2,phase
      complex*16 :: ex
      complex*16, allocatable :: t(:,:),ff(:)
      pi2=pi+pi
      do i=1,nplw_gw
        allocate(t(nr_full,npnt))
        t=(0.d0,0.d0)
        ii=0
        do ind_ista=1,ndim4_nr_red(me_k_b+1)
          ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
          nsta=num_r_star(ista)
          do j0=1,nsta
            ii=ii+1
            j=ind_r_star(j0,ista)
            t(j,:)=pg(i,ii,:)
          enddo
        enddo
        if(nproc_k*nproc_b/=1) call dgop(t,2*nr_full*npnt,'  +',
     &                                   comm_k_b)
        allocate(ff(nr_full))
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          npw=nplwgw(k)
          do j=1,nr_full
            v(:)=r_pnt(1,j)*r0b0(:,1)+r_pnt(2,j)*r0b0(:,2)
     &                               +r_pnt(3,j)*r0b0(:,3)
            phase=pi2*dot_product(pnt(:,k),v)
            ex=dcmplx(cos(phase),sin(phase))
            ff(j)=ex*t(j,k)
          enddo
          call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,ff,1)
          do ind_j=1,nd_b
            j=n0_b+ind_j
            if(j>nplwgw(k)) cycle
            j0=indpw_gw(j,k)
            call zone1_number(gbs(:,j0),rbas,nrdiv,ind)
            pg1(i,ind_j,ind_k)=ff(ind)
          enddo
        enddo
        deallocate(t,ff)
      enddo
      pg1=pg1/nr_full*sqrt(amega)
      end
