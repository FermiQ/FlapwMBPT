      subroutine heg_vrt_for_out(vrt_out,vrt_k_q,vrt_om_nu,vrt_om_q,
     &                           vrt_om_k,vrt_nu_k,vrt_nu_q,ind_q,
     &                           ind_nu,nom,neg,i_wi,lambda_stat_nu,
     &                           lambda_dyn)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: neg
      integer, intent(in) :: ind_q,ind_nu,nom
      complex*16, intent(in) :: lambda_stat_nu(nrr_div_c,nspin),
     &                          lambda_dyn(nrr_div_c,nom,nspin)
      integer, intent(inout) :: i_wi
      complex*16, intent(inout) :: vrt_out(nc_line,0:n_nu),
     &                             vrt_k_q(nc_line,nc_line),
     &                             vrt_nu_q(0:n_nu_exa,nc_line),
     &                        vrt_om_nu(0:n_omega_exa,0:n_nu_exa),
     &                             vrt_om_q(0:n_omega_exa,nc_line),
     &                             vrt_om_k(0:n_omega_exa,nc_line),
     &                             vrt_nu_k(0:n_nu_exa,nc_line)
      integer :: iq,i_nu,i_omega,ig,kf,ik0,ik,ik1,iql
      real*8 :: w1,sp,v(3)
      i_nu=me_t*ndim3_nu+ind_nu-1
      iq=n_mpi_k_red(me_k+1)+ind_q
      iql=k_npnt_in_line_c(iq)
      if(iql/=0) then
c --- Find k=k_f along the same direction as IQ --------------------
        kf=kline_in_npnt_c(n_kl_fermi_c)
        if(iql==1) then
          ik0=kf
        else
          do ig=1,ngroup
            ik=ip_k_c(kf,ig)
            sp=dot_product(q_sph_c(:,ik),q_irr_c(:,iq))
            sp=sp/qabs_star_c(iq)/qabs_star_c(kf)
            if(abs(1.d0-sp)<1.d-5) then
              ik0=ik
              goto 3
            endif
          enddo
        endif
3       ik0=q_sph_in_cube_c(ik0)
c ----------------------- Static Part ------------------------------
        vrt_out(iql,i_nu)=(1.d0,0.d0)+lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
        if(neg) then
          v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
          call zone1_number(v,rb0_c,nr_div_c,ik1)
          vrt_out(iql,i_nu)=vrt_out(iql,i_nu)
     &                     +conjg(lambda_dyn(ik1,i_wi,nspin))
        else
          vrt_out(iql,i_nu)=vrt_out(iql,i_nu)
     &                     +lambda_dyn(ik0,i_wi,nspin)
        endif
c ------------------- VRT_K_Q -------------------------------------
        if(i_nu==0) then
          do ik1=1,nrr_div_c
            ik=q_cube_in_sph_c(ik1)
            if(ik==0) cycle
            ik0=i_kref_c(ik)
            if(iq==1) then
              ik0=k_npnt_in_line_c(ik0)
            else if(ik==1) then
              ik0=1
            else
              sp=dot_product(q_sph_c(:,ik),q_irr_c(:,iq))
              sp=sp/qabs_star_c(iq)/qabs_star_c(ik0)
              if(abs(1.d0-sp)<1.d-7) then
                ik0=k_npnt_in_line_c(ik0)
              else
                ik0=0
              endif
            endif
            if(ik0==0) cycle
c ----------------------- Static Part ------------------------------
            vrt_k_q(ik0,iql)=(1.d0,0.d0)+lambda_stat_nu(ik1,nspin)
c ----------------------- Dynamic Part -----------------------------
            vrt_k_q(ik0,iql)=vrt_k_q(ik0,iql)
     &                      +lambda_dyn(ik1,1,nspin)
          enddo
        endif
c ------------------- VRT_NU_Q -------------------------------------
        if(i_nu<=n_nu_exa) then
          w1=pi/betta_t
          ik0=q_sph_in_cube_c(1)
c ----------------------- Static Part ------------------------------
          vrt_nu_q(i_nu,iql)=(1.d0,0.d0)+lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
          if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
            vrt_nu_q(i_nu,iql)=vrt_nu_q(i_nu,iql)
     &                        +lambda_dyn(ik0,1,nspin)
          else
            v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
            call zone1_number(v,rb0_c,nr_div_c,ik1)
            w1=-w1+w_nu(i_nu)
            w1=w1-w_om_adapt_nu(1,i_nu)
            i_wi=idnint(betta_t*w1/pi/2)+1
            vrt_nu_q(i_nu,iql)=vrt_nu_q(i_nu,iql)
     &                        +conjg(lambda_dyn(ik1,i_wi,nspin))
          endif
        endif
c ------------------- VRT_OM_NU -------------------------------------
        if(iq==1.and.i_nu<=n_nu_exa) then
          ik0=q_sph_in_cube_c(1)
          do i_omega=0,n_omega_exa
            w1=w_omega(i_omega)
c ----------------------- Static Part ------------------------------
            vrt_om_nu(i_omega,i_nu)=(1.d0,0.d0)
     &                             +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_om_nu(i_omega,i_nu)=vrt_om_nu(i_omega,i_nu)
     &                               +lambda_dyn(ik0,i_wi,nspin)
            else
              w1=-w1+w_nu(i_nu)
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_om_nu(i_omega,i_nu)=vrt_om_nu(i_omega,i_nu)
     &                               +conjg(lambda_dyn(ik0,i_wi,nspin))
            endif
          enddo
        endif
c ------------------- VRT_OM_Q -------------------------------------
        if(i_nu==0) then
          ik0=q_sph_in_cube_c(1)
          do i_omega=0,n_omega_exa
c ----------------------- Static Part ------------------------------
            vrt_om_q(i_omega,iql)=(1.d0,0.d0)
     &                           +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            vrt_om_q(i_omega,iql)=vrt_om_q(i_omega,iql)
     &                           +lambda_dyn(ik0,i_omega+1,nspin)
          enddo
        endif
c ------------------- VRT_OM_K -------------------------------------
        if(iq==1.and.i_nu==0) then
          do ik=1,nc_line
            ik0=kline_in_npnt_c(ik)
            ik0=q_sph_in_cube_c(ik0)
            do i_omega=0,n_omega_exa
c ----------------------- Static Part ------------------------------
              vrt_om_k(i_omega,ik)=(1.d0,0.d0)
     &                            +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
              vrt_om_k(i_omega,ik)=vrt_om_k(i_omega,ik)
     &                            +lambda_dyn(ik0,i_omega+1,nspin)
            enddo
          enddo
        endif
c ------------------- VRT_NU_K -------------------------------------
        if(iq==1.and.i_nu<=n_nu_exa) then
          do ik=1,nc_line
            ik0=kline_in_npnt_c(ik)
            ik0=q_sph_in_cube_c(ik0)
c ----------------------- Static Part ------------------------------
            vrt_nu_k(i_nu,ik)=(1.d0,0.d0)
     &                       +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            w1=pi/betta_t
            if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
              vrt_nu_k(i_nu,ik)=vrt_nu_k(i_nu,ik)
     &                         +lambda_dyn(ik0,1,nspin)
            else
              v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
              call zone1_number(v,rb0_c,nr_div_c,ik1)
              w1=-w1+w_nu(i_nu)
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_nu_k(i_nu,ik)=vrt_nu_k(i_nu,ik)
     &                         +conjg(lambda_dyn(ik1,i_wi,nspin))
            endif
          enddo
        endif
      endif
      end
      
      
      
      subroutine heg_vrt_for_out_line(vrt_out,vrt_k_q,vrt_om_nu,
     &                                vrt_om_q,vrt_om_k,vrt_nu_k,
     &                                vrt_nu_q,iql,ind_nu,nom,neg,i_wi,
     &                                lambda_stat_nu,lambda_dyn)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical, intent(in) :: neg
      integer, intent(in) :: iql,ind_nu,nom
      complex*16, intent(in) :: lambda_stat_nu(nrr_div_c,nspin),
     &                          lambda_dyn(nrr_div_c,nom,nspin)
      integer, intent(inout) :: i_wi
      complex*16, intent(inout) :: vrt_out(nc_line,0:n_nu),
     &                             vrt_k_q(nc_line,nc_line),
     &                             vrt_nu_q(0:n_nu_exa,nc_line),
     &                        vrt_om_nu(0:n_omega_exa,0:n_nu_exa),
     &                             vrt_om_q(0:n_omega_exa,nc_line),
     &                             vrt_om_k(0:n_omega_exa,nc_line),
     &                             vrt_nu_k(0:n_nu_exa,nc_line)
      integer :: iq,i_nu,i_omega,ig,kf,ik0,ik,ik1
      real*8 :: w1,sp,v(3)
      i_nu=me_t*ndim3_nu+ind_nu-1
      iq=kline_in_npnt_c(iql)
c --- Find k=k_f along the same direction as IQ --------------------
      kf=kline_in_npnt_c(n_kl_fermi_c)
      if(iql==1) then
        ik0=kf
      else
        do ig=1,ngroup
          ik=ip_k_c(kf,ig)
          sp=dot_product(q_sph_c(:,ik),q_irr_c(:,iq))
          sp=sp/qabs_star_c(iq)/qabs_star_c(kf)
          if(abs(1.d0-sp)<1.d-5) then
            ik0=ik
            goto 3
          endif
        enddo
      endif
3     ik0=q_sph_in_cube_c(ik0)
c ----------------------- Static Part ------------------------------
      vrt_out(iql,i_nu)=(1.d0,0.d0)+lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
      if(neg) then
        v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
        call zone1_number(v,rb0_c,nr_div_c,ik1)
        vrt_out(iql,i_nu)=vrt_out(iql,i_nu)
     &                   +conjg(lambda_dyn(ik1,i_wi,nspin))
      else
        vrt_out(iql,i_nu)=vrt_out(iql,i_nu)
     &                   +lambda_dyn(ik0,i_wi,nspin)
      endif
c ------------------- VRT_K_Q -------------------------------------
      if(i_nu==0) then
        do ik1=1,nrr_div_c
            ik=q_cube_in_sph_c(ik1)
            if(ik==0) cycle
            ik0=i_kref_c(ik)
            if(iq==1) then
              ik0=k_npnt_in_line_c(ik0)
            else if(ik==1) then
              ik0=1
            else
              sp=dot_product(q_sph_c(:,ik),q_irr_c(:,iq))
              sp=sp/qabs_star_c(iq)/qabs_star_c(ik0)
              if(abs(1.d0-sp)<1.d-7) then
                ik0=k_npnt_in_line_c(ik0)
              else
                ik0=0
              endif
            endif
            if(ik0==0) cycle
c ----------------------- Static Part ------------------------------
            vrt_k_q(ik0,iql)=(1.d0,0.d0)+lambda_stat_nu(ik1,nspin)
c ----------------------- Dynamic Part -----------------------------
            vrt_k_q(ik0,iql)=vrt_k_q(ik0,iql)
     &                      +lambda_dyn(ik1,1,nspin)
          enddo
        endif
c ------------------- VRT_NU_Q -------------------------------------
        if(i_nu<=n_nu_exa) then
          w1=pi/betta_t
          ik0=q_sph_in_cube_c(1)
c ----------------------- Static Part ------------------------------
          vrt_nu_q(i_nu,iql)=(1.d0,0.d0)+lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
          if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
            vrt_nu_q(i_nu,iql)=vrt_nu_q(i_nu,iql)
     &                        +lambda_dyn(ik0,1,nspin)
          else
            v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
            call zone1_number(v,rb0_c,nr_div_c,ik1)
            w1=-w1+w_nu(i_nu)
            w1=w1-w_om_adapt_nu(1,i_nu)
            i_wi=idnint(betta_t*w1/pi/2)+1
            vrt_nu_q(i_nu,iql)=vrt_nu_q(i_nu,iql)
     &                        +conjg(lambda_dyn(ik1,i_wi,nspin))
          endif
        endif
c ------------------- VRT_OM_NU -------------------------------------
        if(iq==1.and.i_nu<=n_nu_exa) then
          ik0=q_sph_in_cube_c(1)
          do i_omega=0,n_omega_exa
            w1=w_omega(i_omega)
c ----------------------- Static Part ------------------------------
            vrt_om_nu(i_omega,i_nu)=(1.d0,0.d0)
     &                             +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_om_nu(i_omega,i_nu)=vrt_om_nu(i_omega,i_nu)
     &                               +lambda_dyn(ik0,i_wi,nspin)
            else
              w1=-w1+w_nu(i_nu)
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_om_nu(i_omega,i_nu)=vrt_om_nu(i_omega,i_nu)
     &                               +conjg(lambda_dyn(ik0,i_wi,nspin))
            endif
          enddo
        endif
c ------------------- VRT_OM_Q -------------------------------------
        if(i_nu==0) then
          ik0=q_sph_in_cube_c(1)
          do i_omega=0,n_omega_exa
c ----------------------- Static Part ------------------------------
            vrt_om_q(i_omega,iql)=(1.d0,0.d0)
     &                           +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            vrt_om_q(i_omega,iql)=vrt_om_q(i_omega,iql)
     &                           +lambda_dyn(ik0,i_omega+1,nspin)
          enddo
        endif
c ------------------- VRT_OM_K -------------------------------------
        if(iq==1.and.i_nu==0) then
          do ik=1,nc_line
            ik0=kline_in_npnt_c(ik)
            ik0=q_sph_in_cube_c(ik0)
            do i_omega=0,n_omega_exa
c ----------------------- Static Part ------------------------------
              vrt_om_k(i_omega,ik)=(1.d0,0.d0)
     &                            +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
              vrt_om_k(i_omega,ik)=vrt_om_k(i_omega,ik)
     &                            +lambda_dyn(ik0,i_omega+1,nspin)
            enddo
          enddo
        endif
c ------------------- VRT_NU_K -------------------------------------
        if(iq==1.and.i_nu<=n_nu_exa) then
          do ik=1,nc_line
            ik0=kline_in_npnt_c(ik)
            ik0=q_sph_in_cube_c(ik0)
c ----------------------- Static Part ------------------------------
            vrt_nu_k(i_nu,ik)=(1.d0,0.d0)
     &                       +lambda_stat_nu(ik0,nspin)
c ----------------------- Dynamic Part -----------------------------
            w1=pi/betta_t
            if(w_om_adapt_nu(1,i_nu)-w1<=1.d-12) then
              vrt_nu_k(i_nu,ik)=vrt_nu_k(i_nu,ik)
     &                         +lambda_dyn(ik0,1,nspin)
            else
              v=-q_mesh_heg_c(:,ik0)+q_irr_c(:,iq)
              call zone1_number(v,rb0_c,nr_div_c,ik1)
              w1=-w1+w_nu(i_nu)
              w1=w1-w_om_adapt_nu(1,i_nu)
              i_wi=idnint(betta_t*w1/pi/2)+1
              vrt_nu_k(i_nu,ik)=vrt_nu_k(i_nu,ik)
     &                         +conjg(lambda_dyn(ik1,i_wi,nspin))
            endif
          enddo
        endif
      end