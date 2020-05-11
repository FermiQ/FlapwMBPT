      subroutine polarizability_heg
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical :: w_lda,w_old
      integer :: i_nu,i,iq
      real*8 :: pi2a,hh,q
      real*8, allocatable :: pl(:),spl(:,:)
      pi2a=(pi+pi)/par
! --- Decide if we use W_LDA or not --------------------------------
! --- Decide if we use W_OLD or not --------------------------------
      w_lda=.false.
      w_old=.false.
      if(ubi==' gw') then
        if(w_sc_gw=='non') then
          if(ubi_0=='dft') w_lda=.true.
          if(ubi_0==' gw') w_old=.true.
        endif
      else if(ubi==' qp') then
        if(w_sc_qp=='non') then
          if(ubi_0=='dft') w_lda=.true.
          if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
          if(ubi_0==' qp') w_old=.true.
        endif
      else if(ubi=='psi') then
        if(w_sc_psi=='non') then
          if(ubi_0=='dft') w_lda=.true.
          if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
          if(ubi_0==' qp'.and.w_sc_qp=='non') w_old=.true.
          if(ubi_0=='psi') w_old=.true.
        endif
      else if(ubi=='bsp') then
        if(w_sc_bsp=='non') then
          if(ubi_0=='dft') w_lda=.true.
          if(ubi_0==' gw'.and.w_sc_gw=='non') w_old=.true.
          if(ubi_0==' qp'.and.w_sc_qp=='non') w_old=.true.
          if(ubi_0=='psi'.and.w_sc_psi=='non') w_old=.true.
          if(ubi_0=='bsp') w_old=.true.
        endif
      endif
      if(w_lda.and.bse_kernel_p=='0000') w_lda=.false.
      if(.not.w_old) then
c -------- Polarizability correction ---------------------------
        if(w_lda) then
          call heg_w_dft
          p_q_tau_vertex_heg=0.d0
          p_q_nu_vertex_heg=0.d0
          do i=0,n_q_heg
            call nu_to_tau_spl(p_q_nu_heg(i,:),p_q_tau_heg(i,:),1)
          enddo
          p_q_tau_heg=p_q_tau_heg/2
        else
c -------------- GG - part -------------------------------------------
          call heg1_pol
c --------------------------------------------------------------------          
          if(ubi=='psi'.or.ubi=='bsp') then
            if(ubi=='psi') call heg_vertex_p_psi
            if(ubi=='bsp') call heg_vertex_p_bse
            allocate(pl(nc_line))
            allocate(spl(nc_line-1,4))
            hh=(qabs_line_c(2)-qabs_line_c(1))*pi2a
            do i_nu=0,n_nu
              do iq=1,nc_line
                i=kline_in_npnt_c(iq)
                pl(iq)=p_q_nu_vertex_heg(i_nu,iq)
              enddo
              call spline_hmg(hh,pl,spl(1,1),spl(1,2),spl(1,3),
     &                        spl(1,4),nc_line-1,0,0.d0,0.d0)
              do iq=0,n_q_heg
                q=q_m_heg(iq)
                if(q>hh*(nc_line-1)) then
                  p_q_nu_vertex_rad(iq,i_nu)=0.d0
                else
                  i=q/hh
                  if(i*hh<q) i=i+1
                  if(i==0) i=1
                  q=q-hh*(i-1)
                  p_q_nu_vertex_rad(iq,i_nu)=spl(i,1)+q*(spl(i,2)+q*(
     &                                       spl(i,3)+q*spl(i,4)))
                endif
              enddo
            enddo
            deallocate(pl,spl)
c   ------- Transform NU ---> TAU representation ---------------------
            do i=0,n_q_heg
              call nu_to_tau_spl(p_q_nu_vertex_rad(i,:),
     &                           p_q_tau_vertex_heg(i,:),1)
            enddo
            p_q_tau_vertex_heg=p_q_tau_vertex_heg/2
            call output_heg1_p_nu_vrt
            do i_nu=0,n_nu
              do i=0,n_q_heg
                p_q_nu_heg(i,i_nu)=p_q_nu_heg(i,i_nu)
     &                            +p_q_nu_vertex_rad(i,i_nu)
              enddo
            enddo  !! over i_nu
            do i=0,n_q_heg
              call nu_to_tau_spl(p_q_nu_heg(i,:),p_q_tau_heg(i,:),1)
            enddo
            p_q_tau_heg=p_q_tau_heg/2
          endif
        endif
        if(.not.w_lda) then
          call output_heg1_p_nu(p_q_nu_heg,0)
          call heg1_w
          call output_heg1_p_nu(p_q_nu_heg,3)
        endif
      endif
      end