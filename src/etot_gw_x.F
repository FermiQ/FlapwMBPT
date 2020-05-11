      subroutine etot_gw_x(ispin,k,ind_k,gx_tau,ndim)
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ndim,k,ispin,ind_k
      real*8, intent(in) :: gx_tau(ndim)
      integer :: i,ii,i0,i1
      real*8 :: const,tmp,free_energy_x,s,ee,de,green0
      const=2.d0*wgt(k)/dfloat(nspin*nrel)
      do i=1,ndim
        d_free_x=d_free_x+const
     &                   *free_energy_x(e_bnd(i,k,ispin)-chem_pot)
        if(k_integral=='FD') then
          ebnd_x=ebnd_x-const*gx_tau(i)*e_bnd(i,k,ispin)
        else if(k_integral=='O1') then
          s=0.d0
          do ii=1,nn_k_int
            ee=e_bnd(i,k,ispin)
     &        +de_bnd(i,k,ispin,1)*k_int_vec(1,ii)
     &        +de_bnd(i,k,ispin,2)*k_int_vec(2,ii)
     &        +de_bnd(i,k,ispin,3)*k_int_vec(3,ii)
            de=ee-chem_pot
            s=s+green0(de,betta_t)*ee
          enddo
          s=s/nn_k_int
          ebnd_x=ebnd_x-const*s
        endif
      enddo
      if(ubi/='dft') then
        i0=n0_tb_ij(me_t_b+1)
        i1=i0+nd_tb_ij(me_t_b+1)
        tmp=0.d0
        do i=1,ndim
          ii=i*(i+1)/2
          if(ii<=i0) cycle
          if(ii>i1) cycle
          tmp=tmp-const*gx_tau(i)*sigx_solid_k(ii-i0,ind_k,ispin)
        enddo
        if(nproc_t*nproc_b/=1) call dgop(tmp,1,'  +',comm_t_b)
        gx_sx_val=gx_sx_val+tmp
        tmp=0.d0
        do i=1,ndim
          ii=i*(i+1)/2
          if(ii<=i0) cycle
          if(ii>i1) cycle
          tmp=tmp-const*gx_tau(i)*v_h_bb(ii-i0,ind_k,ispin)
        enddo
        if(nproc_t*nproc_b/=1) call dgop(tmp,1,'  +',comm_t_b)
        gx_vh_val=gx_vh_val+tmp
      endif
      end
