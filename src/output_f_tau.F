      subroutine output_f_tau(f,ch1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      character*6, intent(in) :: ch1
      real*8, intent(in) :: f(nbndf,nbndf,2,ndim3_tau,
     &                        ndim3_k(me_k+1),nspin)
      integer :: k,i_len,n,ispin,ind_k,ind_tau,i_tau,i
      real*8 :: ar,br
      real*8, allocatable :: sig(:,:)
      allocate(sig(0:n_tau,n_bnd_out))
      sig=0.d0
      do n=1,n_bnd_out
        ispin=inf_bnd_out(3,n)
        k=inf_bnd_out(2,n)
        i=n_semic+inf_bnd_out(1,n)
        do ind_k=1,ndim3_k(me_k+1)
          if(k/=n3_mpi_k(me_k+1)+ind_k) cycle
          do ind_tau=1,ndim3_tau
            ar=f(i,i,1,ind_tau,ind_k,ispin)
            br=f(i,i,2,ind_tau,ind_k,ispin)
            i_tau=me_t*ndim3_tau+ind_tau-1
            sig(i_tau,n)=0.5d0*(ar+br)
            sig(n_tau-i_tau,n)=0.5d0*(br-ar)
          enddo
        enddo
      enddo
c --------- Collect from all processes  --------------------------------
      if(nproc_k/=1) then
        call dgop(sig,n_bnd_out*(n_tau+1),'  +',comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(sig,n_bnd_out*(n_tau+1),'  +',comm_t)
      endif
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//ch1//'_band_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//ch1//'_band_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//ch1//'_band_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//ch1//'_band_tau.bsp')
        endif
      endif
      if(maswrk) then
        write(3,'(a9,6(1x,a3,i5))')' Tau     ',
     &                      ('Isp',inf_bnd_out(3,n),n=1,n_bnd_out)
        write(3,'(a9,6(1x,a3,i5))')'         ',
     &                      ('  K',inf_bnd_out(2,n),n=1,n_bnd_out)
        write(3,'(a9,6(1x,a3,i5))')'         ',
     &                      ('BND',inf_bnd_out(1,n),n=1,n_bnd_out)
        do i_tau=0,n_tau
          write(3,'(f9.4,6(1x,f8.4))')tau_mesh(i_tau),
     &                  (sig(i_tau,n),n=1,n_bnd_out)
        enddo
      endif
      if(maswrk) close(3)
      deallocate(sig)
      end
