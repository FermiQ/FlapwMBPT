      subroutine output_sigc_rax(f,ch1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      character*6, intent(in) :: ch1
      real*8, intent(in) :: f(nbndf,2,ndim3_omega,ndim3_k(me_k+1),nspin)
      integer :: k,i_len,n,ispin,ind_k,ind_omega,i_omega,i
      real*8 :: fc0,deltae
      real*8, allocatable :: erax(:)
      complex*16, allocatable :: sig(:,:),sig_out(:,:)
      allocate(sig(0:n_omega,n_bnd_out))
      allocate(sig_out(0:ndos,n_bnd_out))
      allocate(erax(0:ndos))
      deltae=(emaxdos-emindos)/ndos
      do i=0,ndos
        erax(i)=emindos+deltae*i
      enddo
      sig=(0.d0,0.d0)
      do n=1,n_bnd_out
        ispin=inf_bnd_out(3,n)
        k=inf_bnd_out(2,n)
        i=n_semic+inf_bnd_out(1,n)
        do ind_k=1,ndim3_k(me_k+1)
          if(k/=n3_mpi_k(me_k+1)+ind_k) cycle
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
            sig(i_omega,n)=dcmplx(f(i,1,ind_omega,ind_k,ispin),
     &                            f(i,2,ind_omega,ind_k,ispin))
          enddo
        enddo
      enddo
c --------- Collect from all processes  --------------------------------
      if(nproc_k/=1) then
        call dgop(sig,2*n_bnd_out*(n_omega+1),'  +',comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(sig,2*n_bnd_out*(n_omega+1),'  +',comm_t)
      endif
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.gw')
          open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.qp')
          open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.psi')
          open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_Re_'//ch1//'_band_ac.bsp')
          open(4,file=allfile(1:i_len)//'_Im_'//ch1//'_band_ac.bsp')
        endif
      endif
      fc0=evolt/2.d0
      do n=1,n_bnd_out
        call a_cont_ferm(ndos+1,erax,sig(0,n),sig_out(0,n),0.d0)
      enddo
      if(maswrk) then
        write(3,'(a9,6(1x,a3,i5))')' Om (eV) ',
     &                      ('Isp',inf_bnd_out(3,n),n=1,n_bnd_out)
        write(3,'(a9,6(1x,a3,i5))')'         ',
     &                      ('  K',inf_bnd_out(2,n),n=1,n_bnd_out)
        write(3,'(a9,6(1x,a3,i5))')'         ',
     &                      ('BND',inf_bnd_out(1,n),n=1,n_bnd_out)
        write(4,'(a9,6(1x,a3,i5))')' Om (eV) ',
     &                      ('Isp',inf_bnd_out(3,n),n=1,n_bnd_out)
        write(4,'(a9,6(1x,a3,i5))')'         ',
     &                      ('  K',inf_bnd_out(2,n),n=1,n_bnd_out)
        write(4,'(a9,6(1x,a3,i5))')'         ',
     &                      ('BND',inf_bnd_out(1,n),n=1,n_bnd_out)
        do i_omega=0,ndos
          write(3,'(f9.4,6(1x,f8.4))')erax(i_omega)*fc0,
     &                  (real(sig_out(i_omega,n))*fc0,n=1,n_bnd_out)
          write(4,'(f9.4,6(1x,f8.4))')erax(i_omega)*fc0,
     &                  (imag(sig_out(i_omega,n))*fc0,n=1,n_bnd_out)
        enddo
      endif
      if(maswrk) close(3)
      if(maswrk) close(4)
      deallocate(sig_out,sig,erax)
      end
