      subroutine sigma_x_div(ind_k,k,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_k,k,ispin
      integer :: i,j,n,ii
      real*8 :: w,f,ee,fermi_dirac
      complex*16, allocatable :: tmp(:,:)
      if(ubi==' hf'.and.hybrid_type>2) return
      w=wgt(1)*8.d0*pi/amega*q2aver
      n=n_bnd(k,ispin)
      allocate(tmp(n,n))
      tmp=(0.d0,0.d0)
      if(ubi==' hf') then
        do i=1,n
          if(nn_k_int/=1) then
            f=0.d0
            do ii=1,nn_k_int
              ee=e_bnd(i,k,ispin)
     &          +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &          +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &          +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
              f=f+fermi_dirac(ee-chem_pot)
            enddo
            f=-f/nn_k_int
          else
            f=-g_full_00(i,ind_k,ispin)
          endif
          tmp(i,i)=w*f
        enddo
      else
        tmp=(0.d0,0.d0)
        if(me_t==0) then
          call ferm_unpack_tau(tmp,g_full(1,1,1,1,ind_k,ispin),n,n,
     &                         nbndf,2)
          tmp=-tmp
        endif
        if(nproc_t/=1) call dgop(tmp,2*n*n,'  +',comm_t)
        do j=1,n
          do i=1,n
            tmp(i,j)=-w*tmp(i,j)
          enddo
        enddo
      endif
      call pack_hrm_tb(tmp,sigx_solid_k(1,ind_k,ispin),n,n,1)
      deallocate(tmp)
      end
