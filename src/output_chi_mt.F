      subroutine output_chi_mt(ch,chi_mt,chi0_mt,chi_sort,chi0_sort)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*11, intent(in) :: ch
      complex*16, intent(in) :: chi_mt(natom,ndim3_nu),
     &                          chi0_mt(natom,ndim3_nu)
      complex*16, intent(out) :: chi_sort(0:n_nu,nsort),
     &                           chi0_sort(0:n_nu,nsort)
      integer :: i_len,i_nu,i,isort,iatom,ig,jatom,ind_nu
      real*8 :: w,ev
      real*8, allocatable :: tn(:),tt(:)
      chi_sort=(0.d0,0.d0)
      chi0_sort=(0.d0,0.d0)
      do isort=1,nsort
        iatom=iat_1(isort)
        do ig=1,ngroup
          jatom=ip(iatom,ig)
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            chi_sort(i_nu,isort)=chi_sort(i_nu,isort)
     &                          +chi_mt(jatom,ind_nu)
            chi0_sort(i_nu,isort)=chi0_sort(i_nu,isort)
     &                           +chi0_mt(jatom,ind_nu)
          enddo
        enddo
      enddo
      chi_sort=chi_sort/ngroup
      if(nproc_t/=1) call dgop(chi_sort,2*(n_nu+1)*nsort,'  +',
     &                           comm_t)
      chi0_sort=chi0_sort/ngroup
      if(nproc_t/=1) call dgop(chi0_sort,2*(n_nu+1)*nsort,'  +',
     &                           comm_t)
c --------------------------------------------------------------------
      ev=evolt/2.d0
      i_len=len_trim(allfile)
      if(maswrk) then
c --------------------------------------------------------------------
        if(ubi=='dft')open(3,file=allfile(1:i_len)//ch//'Re.dft')
        if(ubi=='dif')open(3,file=allfile(1:i_len)//ch//'Re.dif')
        if(ubi==' hf')open(3,file=allfile(1:i_len)//ch//'Re.hf')
        if(ubi==' qp')open(3,file=allfile(1:i_len)//ch//'Re.qp')
        if(ubi==' gw')open(3,file=allfile(1:i_len)//ch//'Re.gw')
        if(ubi=='psi')open(3,file=allfile(1:i_len)//ch//'Re.psi')
        if(ubi=='bsp')open(3,file=allfile(1:i_len)//ch//'Re.bsp')
        do i_nu=0,n_nu
          w=w_nu(i_nu)*ev
          write(3,100)w,
     &               (real(chi_sort(i_nu,i))/ev,i=1,min(5,nsort))
        enddo
        close(3)
        if(ubi=='dft')open(3,file=allfile(1:i_len)//ch//'Im.dft')
        if(ubi=='dif')open(3,file=allfile(1:i_len)//ch//'Im.dif')
        if(ubi==' hf')open(3,file=allfile(1:i_len)//ch//'Im.hf')
        if(ubi==' qp')open(3,file=allfile(1:i_len)//ch//'Im.qp')
        if(ubi==' gw')open(3,file=allfile(1:i_len)//ch//'Im.gw')
        if(ubi=='psi')open(3,file=allfile(1:i_len)//ch//'Im.psi')
        if(ubi=='bsp')open(3,file=allfile(1:i_len)//ch//'Im.bsp')
        do i_nu=0,n_nu
          w=w_nu(i_nu)*ev
          write(3,100)w,
     &              (imag(chi_sort(i_nu,i))/ev,i=1,min(5,nsort))
        enddo
        close(3)
c --------------------------------------------------------------------
      endif
      allocate(tn(0:n_nu))
      allocate(tt(0:n_tau/2))
      if(ch=='Rsp_ZZ_MT__'.and.maswrk) write(iun,*)' FDT_ZZ:'
      if(ch=='Rsp_+-_MT__'.and.maswrk) write(iun,*)' FDT_XX+YY:'
      do isort=1,nsort
        tn=chi_sort(:,isort)
        call nu_to_tau_aa_r(tn,1,tt,1,1)
        if(maswrk) write(iun,*)txtel(isort),' : ',-tt(0)
      enddo
      if(ch=='Rsp_ZZ_MT__'.and.maswrk) write(iun,*)' FDT_ZZ_0:'
      if(ch=='Rsp_+-_MT__'.and.maswrk) write(iun,*)' FDT_XX+YY_0:'
      do isort=1,nsort
        tn=chi0_sort(:,isort)
        call nu_to_tau_aa_r(tn,1,tt,1,1)
        if(maswrk) write(iun,*)txtel(isort),' : ',-tt(0)
      enddo
      deallocate(tn,tt)
100   format(f11.5,1x,5(e12.5,1x))
      end
