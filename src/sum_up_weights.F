      subroutine sum_up_weights(weigd,e,ispin,ind_k,n)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,ispin,ind_k
      real*8, intent(in) :: e
      real*8, intent(out) :: weigd(n)
      integer :: j,i,k,ii
      real*8 :: gg,f1,f2,d2,ff,ee
      complex*16 :: om0,omi
      weigd=0.d0
      om0=dcmplx(e,0.d0)
      omi=dcmplx(0.d0,e_small)
      d2=e_small**2
      k=n3_mpi_k(me_k+1)+ind_k
      if(ubi=='dft'.or.ubi==' hf') then
        do i=1,n
          f1=(e+chem_pot-e_bnd(i,k,ispin))**2
          f2=1.d0/(f1+d2)
          if(k_integral=='O1') then
            ff=0.d0
            do ii=1,nn_k_int
              ee=e_bnd(i,k,ispin)
     &          +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &          +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &          +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
              f1=(e+chem_pot-ee)**2
              ff=ff+1.d0/(f1+d2)
            enddo
            f2=ff/nn_k_int
          endif
          weigd(i)=wgt(k)*e_small*f2/pi
        enddo
      else if(ubi==' qp') then
        do i=1,n
          gg=imag((1.d0,0.d0)/(om0+omi+chem_pot-e_qp(i,k,ispin)))
          do j=1,n
            weigd(j)=weigd(j)-wgt(k)*gg*abs(q_qp(i,j,ind_k,ispin))**2
     &                       /pi
          enddo
        enddo
      endif
      weigd=weigd*2.d0/dfloat(nspin)/nrel
      end
      
      
      
      subroutine sum_up_weights_x(weigd,e,n0)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n0
      real*8, intent(in) :: e
      real*8, intent(out) :: weigd(n0,npnt,nspin)
      integer :: i,ispin,ind_k,k,n
      real*8 :: gg
      complex*16 :: om
      weigd=0.d0
      om=dcmplx(e,e_small)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          do i=1,n
            gg=imag((1.d0,0.d0)/(om+chem_pot-e_bnd(i,k,ispin)))
            weigd(i,k,ispin)=-wgt(k)*gg/pi
          enddo
        enddo
      enddo
      if(nproc_k/=1) call dgop(weigd,n0*npnt*nspin,'  +',comm_k)
      weigd=weigd*2.d0/dfloat(nspin)/nrel
      end
