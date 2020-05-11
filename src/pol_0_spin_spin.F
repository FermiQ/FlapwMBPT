      subroutine pol_0_spin_spin(isp1,isp2,pmm,pmi,pii)
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: isp1,isp2
      real*8, intent(inout) :: pii(ncmpl*nd_b_pbint(me_b+1),
     &                             ndim3_tn,ndim3_k(me_k+1)),
     &                         pmi(ncmpl*nd_b_pbmt(me_b+1),nplw_gw,
     &                             ndim3_tn,ndim3_k(me_k+1)),
     &                         pmm(ncmpl*ndim_pbmt(me_b+1),ndim3_tn,
     &                             ndim3_k(me_k+1))
      integer :: ind_tau
      do ind_tau=1,ndim3_tau
c	  
c	---------- MT_MT Part ------------------------------------
c	  	    
        call p_mm_sort(pmm,ind_tau,isp1,isp2)
c	  
c	---------- MT_INT Part ------------------------------------
c	  
        call p_im_r(ind_tau,isp1,isp2,pmi)
c	  
c	---------- INT_INT Part ------------------------------------
c	  
        call p_ii_mem1(ind_tau,isp1,isp2,pii)
        if(ind_tau==1)
     &    call timel('** P_II_SOLID finished *************')
      enddo   !! over ind_tau
      end
