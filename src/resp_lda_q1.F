      subroutine resp_lda_q1(ch,n_ix0,ind_ix0)
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
      integer :: n_ch,i_len,i,isort
      real*8 :: h_xi_rax,ev,dsf,dsf0,d,rho
      real*8, allocatable :: wrx(:)
      complex*16, allocatable :: xi0_nu(:,:),xi_nu(:,:),eps1(:),eps10(:)
      ev=evolt/2.d0
      if(ch=='___00') n_ch=1
      if(ch=='___ZZ') n_ch=1
      if(ch=='00_ZZ') n_ch=3
      if(ch=='___XY') n_ch=1
      allocate(xi0_nu(0:nrax_chi,n_ch))
      allocate(xi_nu(0:nrax_chi,n_ch))
      allocate(eps10(0:nrax_chi))
      allocate(eps1(0:nrax_chi))
      xi0_nu=(0.d0,0.d0)
      xi_nu=(0.d0,0.d0)
      eps1=(0.d0,0.d0)
      eps10=(0.d0,0.d0)
      allocate(wrx(0:nrax_chi))
      h_xi_rax=sqrt(freq_chi)/nrax_chi
      do i=0,nrax_chi
        wrx(i)=(i*h_xi_rax)**2
      enddo
c -------------------------------------------------------------------
      call chi_lda_dynamic(ch,wrx,xi0_nu,xi_nu,n_ix0,ind_ix0,n_ch,eps1,
     &                     eps10)
c -------------------------------------------------------------------
      xi_nu=xi_nu*amega
      xi0_nu=xi0_nu*amega
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//ch//'RSP.dft')
        do i=0,nrax_chi
          write(3,100)wrx(i)*ev,xi0_nu(i,1)/ev,xi_nu(i,1)/ev
        enddo
        close(3)
        if(ch=='___00') then
          open(3,file=allfile(1:i_len)//ch//'EPSm1.dft')
          do i=0,nrax_chi
            write(3,100)wrx(i)*ev,eps10(i),eps1(i)
          enddo
          close(3)
          open(3,file=allfile(1:i_len)//ch//'DynStrFact.dft')
          d=dot_product(pnt(:,q_suscept),pnt(:,q_suscept))
          if(iheg) then
            rho=ro_el_gas
          else
            rho=0.d0
            do isort=1,nsort
              rho=rho+z(isort)*nhsort(isort)
            enddo
            rho=rho/amega
          endif
          d=-d/2/par**2/rho
          do i=0,nrax_chi
            dsf0=d*imag(eps10(i))
            dsf=d*imag(eps1(i))
            write(3,100)wrx(i)*ev,dsf0,dsf
          enddo
          close(3)
        endif
      endif
      deallocate(xi0_nu,xi_nu,eps10,eps1,wrx)
      call timel('**** RESP_LDA_Q1 finished **********')
100   format(5(e12.5,1x))
      end
