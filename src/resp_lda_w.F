      subroutine resp_lda_w(ch,n_ch0,ind_ch0,n_ix0,ind_ix0)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*5, intent(in) :: ch
      integer, intent(in) :: n_ch0,ind_ch0(n_ch0),n_ix0,ind_ix0(n_ix0)
      integer :: ind_k
c ------------------------------------------------------------------
      do ind_k=1,ndim3_k(me_k+1)
        if(lsclp) then
#ifdef SCLP
          call p_full_lda_sclp(ch,ind_k,n_ch0,ind_ch0,n_ix0,ind_ix0)
#endif
        else
          call p_full_lda(ch,ind_k,n_ch0,ind_ch0,n_ix0,ind_ix0)
        endif
      enddo  !! over ind_k
      call timel('**** RESP_LDA_W finished ***********')
      end
