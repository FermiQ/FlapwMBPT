      subroutine pol_sclp_spin_spin(isp1,isp2,pwsclp)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: ind_tau
      do ind_tau=1,ndim3_tau
c
c   ---------- MT_MT Part ------------------------------------
c
        call p_mm_sclp(pwsclp,ind_tau,isp1,isp2)
c
c   ---------- MT_INT Part ------------------------------------
c
        call p_im_sclp(ind_tau,isp1,isp2,pwsclp)
c
c   ---------- INT_INT Part ------------------------------------
c
        call p_ii_sclp(ind_tau,isp1,isp2,pwsclp)
        if(ind_tau==1)
     &    call timel('** P_II_SOLID finished *************')
      enddo   !! over ind_tau
      end
