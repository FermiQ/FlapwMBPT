      subroutine optics_gen_no_lcf(nrax,wrx,e)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use units_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nrax
      complex*16, intent(in) :: wrx(0:nrax)
      complex*16, intent(out) :: e(3,3,0:nrax)
      integer :: k,ispin,n,i,j,ind_k,ig,i_nu,ind_nu,ind_tau,isort,
     &           i_tau
      real*8 :: pi8,tu
      complex*16 :: pp(3,3),p1(3,3),p2(3,3),per(3,3),cc,rho0,wp2
      real*8, allocatable :: wr(:),gx(:,:)
      complex*16, allocatable :: ea(:,:,:),p_tau(:,:,:,:),
     &                           tmp(:,:,:),aa(:,:,:),bb(:,:,:),
     &                           p_nu(:,:,:),ff(:),xir(:,:,:)
      pi8=8.d0*pi
      allocate(p_tau(3,3,2,ndim3_tau))
      p_tau=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(ea(n,n,3))
          call bands_gradient(nbndf,n,k,ispin,z_bnd(1,1,ind_k,ispin),
     &                        ev_bnd(1,1,ind_k,ispin),ea)
          allocate(tmp(n,n,2))
          allocate(aa(n,n,3))
          allocate(bb(n,n,3))
          allocate(gx(n,2))
          do ind_tau=1,ndim3_tau
            if(ubi=='dft'.or.ubi==' hf') then
              i_tau=me_t*ndim3_tau+ind_tau-1
              tu=tau_mesh(i_tau)
              call g_x_tau(ispin,k,gx(1,1),tu,n,chem_pot)
              tu=betta_t-tu
              call g_x_tau(ispin,k,gx(1,2),tu,n,chem_pot)
              do j=1,3
                do i=1,n
                  aa(i,:,j)=gx(i,2)*ea(i,:,j)
                  bb(:,i,j)=gx(i,1)*ea(i,:,j)
                enddo
              enddo
            else
              call ferm_unpack_tau(tmp,
     &                             g_full(1,1,1,ind_tau,ind_k,ispin),n,
     &                             n,nbndf,3)
              do i=1,3
                call zgemm('n','n',n,n,n,(1.d0,0.d0),tmp(1,1,2),n,
     &                     ea(1,1,i),n,(0.d0,0.d0),aa(1,1,i),n)
                call zgemm('t','t',n,n,n,(1.d0,0.d0),ea(1,1,i),n,
     &                     tmp(1,1,1),n,(0.d0,0.d0),bb(1,1,i),n)
              enddo
            endif
            call zgemm('t','n',3,3,n*n,(1.d0,0.d0),aa,n*n,bb,n*n,
     &                 (0.d0,0.d0),pp,3)
            p_tau(:,:,1,ind_tau)=p_tau(:,:,1,ind_tau)+wgt(k)*pp
          enddo   !! over ind_tau
          deallocate(ea,tmp,aa,bb,gx)
        enddo  !! over ind_k
      enddo    !! over ispin
      if(nproc_k/=1) call DGOP(p_tau,36*ndim3_tau,'  +',comm_k)
      if(nspin==1.and.irel/=2) p_tau=2.d0*p_tau
      p_tau=4.d0*p_tau/amega
c ---------------- Symmetrization ------------------------------
      do ind_tau=1,ndim3_tau
        per=(0.d0,0.d0)
        do ig=1,ngroup
          do j=1,3
            call rotate_c(p_tau(j,:,1,ind_tau),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_c(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        p_tau(:,:,1,ind_tau)=per
      enddo
c -------- Get P(beta-tau) -------------------------------------
      do ind_tau=1,ndim3_tau
        p_tau(:,:,2,ind_tau)=transpose(p_tau(:,:,1,ind_tau))
      enddo
c --------------------------------------------------------------
      allocate(p_nu(3,3,ndim3_nu))
      call tau_to_nu_gen_spl(p_nu,p_tau,9)
      deallocate(p_tau)
c --------------------------------------------------------------
      allocate(xir(3,3,0:nrax))
      allocate(wr(0:nrax))
      allocate(ff(0:n_nu))
      wr=real(wrx)
      do j=1,3
        do i=1,3
          ff=(0.d0,0.d0)
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            ff(i_nu)=p_nu(i,j,ind_nu)
          enddo
          if(nproc_t/=1) call dgop(ff,2*(n_nu+1),'  +',comm_t)
          tu=maxval(abs(ff))
          if(tu<1.d-12) then
            xir(i,j,:)=(0.d0,0.d0)
          else
            call pade_boson_all(ff,wr,nrax,xir(i,j,:),2)
          endif
        enddo
      enddo
      deallocate(ff,p_nu)
c ----------------------------------------------------------------
      rho0=0.d0
      do isort=1,nsort
        rho0=rho0+z(isort)*nhsort(isort)
      enddo
      rho0=rho0/amega
      wp2=16.d0*pi*rho0
      e=(0.d0,0.d0)
      do i_nu=0,nrax
        cc=(1.d0,0.d0)
        if(i_nu==0) cycle
c        cc=cc-wp2/wr(i_nu)**2
        do i=1,3
          e(i,i,i_nu)=cc
        enddo
        cc=pi8/wr(i_nu)**2
        do j=1,3
          do i=1,3
            e(i,j,i_nu)=e(i,j,i_nu)-cc*xir(i,j,i_nu)
          enddo
        enddo
      enddo
      deallocate(xir,wr)
      end
