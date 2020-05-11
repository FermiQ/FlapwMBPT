      subroutine heg_wi_prepare(i_nu,nom,iq,lambda_dyn,lambda_stat)
c ------------------ Full WI checking ------------------- 
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
      integer, intent(in) :: i_nu,iq,nom
      complex*16, intent(in) :: lambda_stat(nrr_div_c),
     &                          lambda_dyn(nrr_div_c,nom)
      integer :: i_omega,i,k,k0
      complex*16 :: cc
c --------- Density screened vertex ------------------------------      
      do k=1,nrr_div_c
        k0=k
        do i_omega=1,nom
          cc=lambda_stat(k0)+lambda_dyn(k0,i_omega)
          vrt3s_heg(k,i_omega)=(1.d0+cc)/edd_q_nu_heg(i_nu,iq)
     &	                    -(1.d0,0.d0)
        enddo
      enddo
c --------- Density-Current Dielectrix matrix -----------------  
      edc_q_nu_heg(iq,:,i_nu)=-v_q_heg(iq)*pdc_q_nu_heg(iq,:,i_nu)
c --------- Current screened vertex ------------------------------      
      do k=1,nrr_div_c
        do i_omega=1,nom
          do i=1,3
            vrt3s_cur_heg(i,k,i_omega)=vrt3s_cur_heg(i,k,i_omega)
     &	      -edc_q_nu_heg(iq,i,i_nu)*(0.d0,1.d0)*
     &        *(1.d0+vrt3s_heg(k,i_omega))/edd_q_nu_heg(i_nu,iq)
          enddo
        enddo
      enddo
      end