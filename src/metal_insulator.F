      subroutine metal_insulator(ph)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: ph(n_k_head_0,ndim3_nu)
      integer :: ind_nu,i_nu
      real*8 :: p0(0:1)
c ---------- Metal or Insulator ------------------------------------
      p0=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        if(i_nu==0) p0(0)=ph(1,ind_nu)
        if(i_nu==1) p0(1)=ph(1,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(p0,2,'  +',comm_t)
      p0=p0+p_q0_vrt
      metal=.false.
      if(abs(p0(0)-p0(1))>abs(p0(1))) metal=.true.
      if(maswrk) then
        if(metal) write(iun,*)
     &     ' The solid has been identified as METAL'
        if(.not.metal) write(iun,*)
     &     ' The solid has been identified as INSULATOR'
      endif
      end
