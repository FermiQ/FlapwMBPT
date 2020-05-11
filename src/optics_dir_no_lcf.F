      subroutine optics_dir_no_lcf(nrax,wrx,e)
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
     &           i_tau,ii,jj
      real*8 :: pi8,tu,de,fermi_dirac,dd,dd0,dd1,ddr,ddi
      complex*16 :: pp(3,3),p1(3,3),p2(3,3),per(3,3),cc,rho0,wp2
      real*8, allocatable :: wr(:),gx(:,:),fk(:)
      complex*16, allocatable :: ea(:,:,:),
     &                           tmp(:,:,:),aa(:,:,:),bb(:,:,:),
     &                           p_nu(:,:,:),ff(:),xir(:,:,:)
      pi8=8.d0*pi
      allocate(wr(0:nrax))
      wr=real(wrx)
      allocate(xir(3,3,0:nrax))
      xir=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(ea(n,n,3))
          call bands_gradient(nbndf,n,k,ispin,z_bnd(1,1,ind_k,ispin),
     &                        ev_bnd(1,1,ind_k,ispin),ea)
          allocate(fk(n))
          do i=1,n
            de=e_bnd(i,k,ispin)-chem_pot
            fk(i)=fermi_dirac(de)
          enddo
          do i_nu=0,nrax
            pp=(0.d0,0.d0)
            do j=1,n
              do i=1,n



                dd0=e_bnd(j,k,ispin)-e_bnd(i,k,ispin)
                if(abs(dd0)<1.d-8) cycle
                dd=dd0+wr(i_nu)
                dd1=dd*dd+e_small_bos**2
                ddr=dd/dd1
                ddi=-e_small_bos/dd1/dd0**2



                de=fk(j)-fk(i)
                cc=de*dcmplx(ddr,ddi)
c                cc=e_bnd(j,k,ispin)-e_bnd(i,k,ispin)+wrx(i_nu)
c                cc=de/cc
                do jj=1,3
                  do ii=1,3
                    pp(ii,jj)=pp(ii,jj)+cc*ea(j,i,ii)*ea(i,j,jj)
                  enddo
                enddo
              enddo
            enddo
            xir(:,:,i_nu)=xir(:,:,i_nu)+wgt(k)*pp
          enddo   !! over i_nu
          deallocate(ea,fk)
        enddo  !! over ind_k
      enddo    !! over ispin
      if(nproc_k/=1) call DGOP(xir,18*(nrax+1),'  +',comm_k)
      if(nspin==1.and.irel/=2) xir=2.d0*xir
      xir=-4.d0*xir/amega
c ---------------- Symmetrization ------------------------------
      do i_nu=0,nrax
        per=(0.d0,0.d0)
        do ig=1,ngroup
          do j=1,3
            call rotate_c(xir(j,:,i_nu),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_c(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        xir(:,:,i_nu)=per
      enddo
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
        cc=pi8 !/wr(i_nu)**2
        do j=1,3
          do i=1,3
            e(i,j,i_nu)=e(i,j,i_nu)-cc*xir(i,j,i_nu)
          enddo
        enddo
      enddo
      deallocate(xir,wr)
      end
