      subroutine tos_heg1(x,ttos,key)
c     key=0 for free electrons/DFT
c     key=1 for exchange only
c     key=2 for GW, PSI, or BSP
c     key=3 for QP
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(in) :: x
      real*8, intent(inout) :: ttos
      integer :: i,i_omega,ispin
      real*8 :: fermi_dirac,de,ffd,k2,k
      complex*16 :: gxx,gxs
      real*8, allocatable :: tmp(:)
      complex*16, allocatable :: gc(:)
      if(key==0) then
c ---------- DFT -------------------------------
        if(dim_heg==2) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k=k_m_heg(i)
              de=k*k+vxc_heg(ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        else if(dim_heg==3) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k2=k_m_heg(i)**2
              de=k2+vxc_heg(ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k2*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        endif
      else if(key==1) then
c ---------- GW - Exchange ------------------------------
        if(dim_heg==2) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k=k_m_heg(i)
              de=e_star_x(i,ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        else if(dim_heg==3) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k2=k_m_heg(i)**2
              de=e_star_x(i,ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k2*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        endif
      else if(key==2) then
c ---------- GW - Exchange ------------------------------
        if(dim_heg==2) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k=k_m_heg(i)
              de=e_star_x(i,ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        else if(dim_heg==3) then
          do ispin=1,nspin
            do i=0,n_k_heg
              k2=k_m_heg(i)**2
              de=e_star_x(i,ispin)-x
              ffd=fermi_dirac(de)
              ttos=ttos+k2*ffd*k_w_heg(i)
            enddo  !! over i
          enddo
        endif
c ---------- GW - Correlation ------------------------------
        allocate(gc(0:n_omega))
        allocate(tmp(0:n_tau))
        do ispin=1,nspin
          do i=0,n_k_heg
            k2=k_m_heg(i)**2
            do i_omega=0,n_omega
              de=e_star_x(i,ispin)-x
              gxx=(1.d0,0.d0)/dcmplx(-de,w_omega(i_omega))
              gxs=gxx*sigc_k_omega_heg(i,i_omega,ispin)
              gc(i_omega)=gxs*gxx/((1.d0,0.d0)-gxs)
            enddo  !! over i_omega
            tmp(n_tau)=0.5d0*(
     &            dot_product(imag(gc),tau_from_omega(:,0,2))
     &           -dot_product(real(gc),tau_from_omega(:,0,1)))
            ttos=ttos-k2*tmp(n_tau)*k_w_heg(i)
          enddo  !! over i
        enddo 
        deallocate(gc,tmp)
      else if(key==3) then
        call bands_qp_heg1(x)
        do ispin=1,nspin
          do i=0,n_k_heg
            k2=k_m_heg(i)**2
            de=e_star_qp(i,ispin)-x
            ffd=fermi_dirac(de)
            ttos=ttos+k2*ffd*k_w_heg(i)
          enddo  !! over i
        enddo
      endif
      ttos=ttos/2.d0/pi/pi/ro_el_gas
      if(nspin==1) ttos=2.d0*ttos   !! spin
      end