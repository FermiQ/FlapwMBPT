#ifdef SCLP
      subroutine sclp_pw_fin
      use scalapack_mod
      implicit none
      deallocate(ncol_pw,nrow_pw,refrow_pw,refcol_pw,desc_pw)
      call blacs_gridexit(sctxt_pw)
      end
#endif
