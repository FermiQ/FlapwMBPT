      subroutine heg_splines
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
      integer :: i_tau,i_omega,i,ispin,iq,ik
      real*8 :: de,green0,pi2a,pi2,ak,x,dn,sx,k2,k,gc
      complex*16 :: sc,gg
      real*8, allocatable :: st(:)
      complex*16, allocatable :: sk(:)
      pi2=pi+pi
      pi2a=pi2/par
c -------- Splines for SIG_X (k) --------------------------------
      dn=8.d0*pi*ro_el_gas/k_m_heg(n_k_heg)**3
      if(abs(sigx_k_heg(0,1))<1.d-12) dn=0.d0
      do ispin=1,nspin
        call spline_inhmg(k_m_heg,sigx_k_heg(0,ispin),
     &                    sigx_k_spl(1,1,ispin),
     &                    sigx_k_spl(1,2,ispin),
     &                    sigx_k_spl(1,3,ispin),
     &                    sigx_k_spl(1,4,ispin),n_k_heg,2,0.d0,dn)
      enddo
c -------- Splines for SIG (k,omega) --------------------------------
      allocate(sk(0:n_k_heg))
      do ispin=1,nspin
        do i_omega=0,n_omega
          sk=sigx_k_heg(:,ispin)+sigc_k_omega_heg(:,i_omega,ispin)
          call spline_inhmg_c(k_m_heg,sk,
     &                        sig_k_omega_spl(1,1,i_omega,ispin),
     &                        sig_k_omega_spl(1,2,i_omega,ispin),
     &                        sig_k_omega_spl(1,3,i_omega,ispin),
     &                        sig_k_omega_spl(1,4,i_omega,ispin),
     &                        n_k_heg,0,(0.d0,0.d0),(0.d0,0.d0))
          sk=sigc_k_omega_heg(:,i_omega,ispin)
          if(ubi=='psi'.or.ubi=='bsp') sk=
     &                     sk-sig_k_omega_vertex_heg(:,i_omega,ispin)
          call spline_inhmg_c(k_m_heg,sk,
     &                        sigc_gw_k_omega_spl(1,1,i_omega,ispin),
     &                        sigc_gw_k_omega_spl(1,2,i_omega,ispin),
     &                        sigc_gw_k_omega_spl(1,3,i_omega,ispin),
     &                        sigc_gw_k_omega_spl(1,4,i_omega,ispin),
     &                        n_k_heg,0,(0.d0,0.d0),(0.d0,0.d0))
        enddo
      enddo
      deallocate(sk)
c -------- Splines for GC(k,tau) --------------------------------        
      do ispin=1,nspin
        do i_tau=0,n_tau
          call spline_inhmg(k_m_heg,gc_k_tau_heg(0,i_tau,ispin),
     &                      gc_k_tau_spl(1,1,i_tau,ispin),
     &                      gc_k_tau_spl(1,2,i_tau,ispin),
     &                      gc_k_tau_spl(1,3,i_tau,ispin),
     &                      gc_k_tau_spl(1,4,i_tau,ispin),n_k_heg,0,
     &                      0.d0,0.d0)
        enddo
      enddo
c -------- Splines for G(k,tau) --------------------------------        
      do ispin=1,nspin
        do i_tau=0,n_tau
          call spline_inhmg(k_m_heg,g_k_tau_heg(0,i_tau,ispin),
     &                      g_k_tau_spl(1,1,i_tau,ispin),
     &                      g_k_tau_spl(1,2,i_tau,ispin),
     &                      g_k_tau_spl(1,3,i_tau,ispin),
     &                      g_k_tau_spl(1,4,i_tau,ispin),n_k_heg,0,
     &                      0.d0,0.d0)
        enddo
        do i_tau=1,3
          call spline_inhmg(k_m_heg,dg_k_tau_heg(0,i_tau,ispin),
     &                      dg_k_tau_spl(1,1,i_tau,ispin),
     &                      dg_k_tau_spl(1,2,i_tau,ispin),
     &                      dg_k_tau_spl(1,3,i_tau,ispin),
     &                      dg_k_tau_spl(1,4,i_tau,ispin),n_k_heg,0,
     &                      0.d0,0.d0)
        enddo
      enddo
c -------- Splines for G(k,tau)*SIG(k,tau) ------------------------
      allocate(sk(0:n_omega))
      allocate(st(0:n_tau))
      st=0.d0
      do ispin=1,nspin
        do ik=0,n_k_heg
          k=k_m_heg(ik)
          k2=k*k
          sx=sigx_k_heg(ik,ispin)
          de=k2+sx-chem_pot
          do i_omega=0,n_omega
            sc=sigc_k_omega_heg(ik,i_omega,ispin)
            gg=(1.d0,0.d0)/(dcmplx(-de,w_omega(i_omega))-sc)
            sk(i_omega)=sc*gg
            gsig_k_omega_heg(ik,i_omega,ispin)=(sx+sc)*gg
          enddo
c          if(ubi/='dft'.and.ubi/=' hf') call from_omega_to_tau_v(1,sk,
c     &                                                           st)
          do i_tau=0,n_tau
            gsig_k_tau_heg(ik,i_tau,ispin)=st(i_tau)
     &                            +sx*g_k_tau_heg(ik,i_tau,ispin)
          enddo
        enddo
        do i_tau=0,n_tau
          call spline_inhmg(k_m_heg,gsig_k_tau_heg(0,i_tau,ispin),
     &                      gsig_k_tau_spl(1,1,i_tau,ispin),
     &                      gsig_k_tau_spl(1,2,i_tau,ispin),
     &                      gsig_k_tau_spl(1,3,i_tau,ispin),
     &                      gsig_k_tau_spl(1,4,i_tau,ispin),n_k_heg,0,
     &                      0.d0,0.d0)
        enddo
      enddo
      deallocate(sk,st)
! -------------- Asymptotic K --> Inf for SIGC -----------------
      do ispin=1,nspin
        do i_omega=0,n_omega
          sigc_k_omega_as(i_omega,ispin)=
     &      sigc_k_omega_heg(n_k_heg,i_omega,ispin)*k_infinum**3
        enddo
      enddo
! ------- Quantities in CUBE for Vertex part -------------------      
      if(iter_psi+iter_bsp/=0) then
        do ispin=1,nspin
          do iq=1,nstar_c
            k=pi2a*qabs_star_c(iq)
            k2=k*k
            x=x_fermi+asinh((k-k_fermi)/alf_k_heg)
            i=x/step_k_x
            if(i*step_k_x<x) i=i+1
            if(i==0) i=1
            ak=k-k_m_heg(i-1)
            sx=sigx_k_spl(i,1,ispin)+ak*(sigx_k_spl(i,2,ispin)+ak*(
     &         sigx_k_spl(i,3,ispin)+ak*sigx_k_spl(i,4,ispin)))
            e_star_xc(iq,ispin)=k2+sx
            sigx_k_heg_c(iq,ispin)=sx
! -------------------------------------------------------------            
            de=e_star_xc(iq,ispin)-chem_pot
            do i_tau=0,n_tau
              gc=gc_k_tau_spl(i,1,i_tau,ispin)+ak*(
     &           gc_k_tau_spl(i,2,i_tau,ispin)+ak*(
     &           gc_k_tau_spl(i,3,i_tau,ispin)+ak*
     &           gc_k_tau_spl(i,4,i_tau,ispin)))
              g_k_tau_heg_c(i_tau,iq,ispin)=
     &          green0(de,tau_mesh(i_tau))+gc
            enddo
            do i_tau=0,n_tau
              gx_k_tau_heg_c(i_tau,iq,ispin)=
     &           green0(de,tau_mesh(i_tau))
            enddo
            do i_tau=1,3
              gc=dg_k_tau_spl(i,1,i_tau,ispin)+ak*(
     &           dg_k_tau_spl(i,2,i_tau,ispin)+ak*(
     &           dg_k_tau_spl(i,3,i_tau,ispin)+ak*
     &           dg_k_tau_spl(i,4,i_tau,ispin)))
              dg_k_tau_heg_c(i_tau,iq,ispin)=gc
            enddo
! -------------------------------------------------------------
            do i_omega=0,n_omega
              gx_k_omega_heg_c(i_omega,iq,ispin)=(1.d0,0.d0)
     &                             /dcmplx(-de,w_omega(i_omega))
              sc=sig_k_omega_spl(i,1,i_omega,ispin)+ak*(
     &           sig_k_omega_spl(i,2,i_omega,ispin)+ak*(
     &           sig_k_omega_spl(i,3,i_omega,ispin)+ak*
     &           sig_k_omega_spl(i,4,i_omega,ispin)))
              g_k_omega_heg_c(i_omega,iq,ispin)=(1.d0,0.d0)
     &            /(dcmplx(chem_pot-k2,w_omega(i_omega))-sc)
              sigc_k_omega_heg_c(i_omega,iq,ispin)=sc-sx
              sc=sigc_gw_k_omega_spl(i,1,i_omega,ispin)+ak*(
     &           sigc_gw_k_omega_spl(i,2,i_omega,ispin)+ak*(
     &           sigc_gw_k_omega_spl(i,3,i_omega,ispin)+ak*
     &           sigc_gw_k_omega_spl(i,4,i_omega,ispin)))
              sigc_k_omega_heg_c(i_omega,iq,ispin)=sc
            enddo
          enddo
        enddo
      endif
      end