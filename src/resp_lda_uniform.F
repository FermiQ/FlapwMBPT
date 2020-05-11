      subroutine resp_lda_uniform(ch,n_ix0,ind_ix0)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ix0,ind_ix0(n_ix0)
      integer :: n_ch
      real*8 :: xi
      complex*16 :: eps1,eps10
      real*8, allocatable :: xi0_nu(:),xi_nu(:)
      if(ch=='___00') n_ch=1
      if(ch=='___ZZ') n_ch=1
      if(ch=='00_ZZ') n_ch=3
      if(ch=='___XY') n_ch=1
      allocate(xi0_nu(n_ch))
      allocate(xi_nu(n_ch))
      xi0_nu=0.d0
      xi_nu=0.d0
      eps1=0.d0
      eps10=0.d0
c -------------------------------------------------------------------
      call chi_lda_uniform(ch,xi0_nu,xi_nu,n_ix0,ind_ix0,n_ch,eps1,
     &                     eps10)
c -------------------------------------------------------------------
      if(nproc_k/=1) then
        call dgop(xi0_nu,n_ch,'  +',comm_k)
        call dgop(xi_nu,n_ch,'  +',comm_k)
        call dgop(eps10,1,'  +',comm_k)
        call dgop(eps1,1,'  +',comm_k)
      endif
      xi_nu=xi_nu*amega
      xi0_nu=xi0_nu*amega
      if(maswrk) then
        if(ch=='___00') then
          write(iun,*)'Chi_0 Chi Uniform:',xi0_nu(1),xi_nu(1)
        else if(ch=='___ZZ') then
          write(iun,*)' Uniform ZZ Chi_0:'
          write(iun,*)' Volume: ',-3.34588d-4*xi0_nu(1)/amega
          xi=-2.98575d-5*xi0_nu(1)
          write(iun,*)'   Mass: ',xi/molar_mass,' cm^3/g'
          write(iun,*)'  Molar: ',xi,' cm^3/mol'
          write(iun,*)' Uniform ZZ Chi:'
          write(iun,*)' Volume: ',-3.34588d-4*xi_nu(1)/amega
          xi=-2.98575d-5*xi_nu(1)
          write(iun,*)'   Mass: ',xi/molar_mass,' cm^3/g'
          write(iun,*)'  Molar: ',xi,' cm^3/mol'
        endif
	  endif
      deallocate(xi0_nu,xi_nu)
      call timel('**** RESP_LDA_UNIFORM finished *****')
      end
