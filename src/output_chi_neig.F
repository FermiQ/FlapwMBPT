      subroutine output_chi_neig(ch,chi_neig)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*11, intent(in) :: ch
      complex*16, intent(in) :: chi_neig(max_sites,nsort,ndim3_nu)
      integer :: i_nu,isort,iatom,ind_nu,nei
      real*8 :: ev
      real*8, allocatable :: tn(:),tt(:)
      complex*16, allocatable :: chi(:,:)
c --------------------------------------------------------------------
      ev=evolt/2.d0
      allocate(tn(0:n_nu))
      allocate(tt(0:n_tau/2))
      if(ch=='Rsp_ZZ_MT__'.and.maswrk)
     &  write(iun,*)' ZZ: Isort  Neighbor   Nu=0    Tau=0'
      if(ch=='Rsp_+-_MT__'.and.maswrk)
     &  write(iun,*)' XX+YY: Isort  Neighbor   Nu=0    Tau=0'
      do isort=1,nsort
        iatom=iat_1(isort)
        do nei=1,nn_sites(iatom)
          tn=0.d0
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            tn(i_nu)=chi_neig(nei,isort,ind_nu)
          enddo
          if(nproc_t/=1) call dgop(tn,n_nu+1,'  +',comm_t)
          call nu_to_tau_aa_r(tn,1,tt,1,1)
          if(maswrk) write(iun,*)txtel(isort),nei,-tn(0)/ev,-tt(0)
        enddo
      enddo
      deallocate(tn,tt)
      allocate(chi(0:n_nu,nsort))
      chi=(0.d0,0.d0)
      do isort=1,nsort
        do ind_nu=1,ndim3_nu
          i_nu=me_t*ndim3_nu+ind_nu-1
          chi(i_nu,isort)=chi_neig(1,isort,ind_nu)
        enddo
      enddo
      if(nproc_t/=1) call dgop(chi,2*(n_nu+1)*nsort,'  +',comm_t)
      call from_nu_to_rax_mt('ZN',chi)
      deallocate(chi)
      end
