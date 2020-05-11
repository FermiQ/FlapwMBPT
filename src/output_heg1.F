      subroutine output_heg1_g_tau
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
      integer :: i_len,i_tau,i
      character*3 :: ch
      real*8, allocatable :: s(:,:),sc(:,:)
      allocate(s(0:n_tau,0:n_k_heg))
      do i=0,n_k_heg
        s(:,i)=g_k_tau_heg(i,:,nspin)
      enddo
      allocate(sc(0:n_tau,0:n_k_heg))
      do i=0,n_k_heg
        sc(:,i)=gc_k_tau_heg(i,:,nspin)
      enddo
      if(ubi=='dft') ch='dft'
      if(ubi==' hf') ch='hf '
      if(ubi==' qp') ch='qp '
      if(ubi==' gw') ch='gw '
      if(ubi=='psi') ch='psi'
      if(ubi=='bsp') ch='bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_G_tau_k.'//ch)
        do i_tau=0,n_tau
          write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &      	                        (s(i_tau,i),i=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_G_k_tau.'//ch)
        do i=0,n_k_heg
          write(3,'(f9.4,5(1x,f12.7))')k_m_heg(i),
     &	                             (s(i_tau,i),i_tau=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_GC_tau_k.'//ch)
        do i_tau=0,n_tau
          write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &	  	                        (sc(i_tau,i),i=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_GC_k_tau.'//ch)
        do i=0,n_k_heg
          write(3,'(f9.4,5(1x,f12.7))')k_m_heg(i),
     &	                             (sc(i_tau,i),i_tau=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_G_r_tau.'//ch)
        do i=0,n_r_heg
          write(3,'(f9.3,5(1x,f12.7))')r_m_heg(i),
     &	 	             (g_r_tau_heg(i,i_tau,nspin),i_tau=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_G_tau_r.'//ch)
        do i_tau=0,n_tau
          write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &		                   (g_r_tau_heg(i,i_tau,nspin),i=0,4)
        enddo
        close(3)
      endif
      deallocate(s,sc)
      end
      
      
      
      subroutine output_heg1_sig_c_tau
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
      integer :: i_len,i_tau,i,ik
      if(maswrk) then
        i_len=len_trim(allfile)
c ------------ Small K ---------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k.bsp')
        endif
        do i_tau=0,n_tau
          write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &		(sigc_k_tau_heg(i,i_tau,nspin)*evolt/2,i=0,4)
        enddo
        close(3)
c ------------ Large K ---------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k_large.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k_large.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k_large.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_tau_k_large.bsp')
        endif
        do i_tau=0,n_tau
          write(3,'(f9.4,4(1x,f13.7))')tau_mesh(i_tau),
     &    (sigc_k_tau_heg(i,i_tau,nspin)*evolt/2,i=n_k_heg-3,n_k_heg)
        enddo
        close(3)
c -------------------------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_C_K_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_K_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_C_K_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_C_K_tau.bsp')
        endif
        do ik=0,n_k_heg
          write(3,'(f9.3,5(1x,f12.7))')k_m_heg(ik),
     &		(sigc_k_tau_heg(ik,i,nspin)*evolt/2,i=0,4)
        enddo
        close(3)
      endif
      end 
      
      
      
      subroutine output_heg1_sig_omega
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
      integer :: i_len,i_omega,i
      complex*16, allocatable :: s(:,:),sv(:,:)
      allocate(s(0:n_omega,0:n_k_heg))
      allocate(sv(0:n_omega,0:n_k_heg))
      do i=0,n_k_heg
        s(:,i)=sigc_k_omega_heg(i,:,nspin)
      enddo
      if(ubi=='psi'.or.ubi=='bsp') then
        do i=0,n_k_heg
          sv(:,i)=sig_k_omega_vertex_heg(i,:,nspin)
        enddo
      endif
      if(maswrk) then
        i_len=len_trim(allfile)
c ------------ Small K -----------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re.psi')
          open(4,file=allfile(1:i_len)//'_SIGC_VRT_omega_Re.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re.bsp')
          open(4,file=allfile(1:i_len)//'_SIGC_VRT_omega_Re.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=0,4)
        enddo
        if(ubi=='psi'.or.ubi=='bsp') then
          do i_omega=0,n_omega
            write(4,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &          (real(sv(i_omega,i))*evolt/2,i=0,4)
          enddo
        endif
        close(3)
        if(ubi=='psi'.or.ubi=='bsp') close(4)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im.psi')
          open(4,file=allfile(1:i_len)//'_SIGC_VRT_omega_Im.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im.bsp')
          open(4,file=allfile(1:i_len)//'_SIGC_VRT_omega_Im.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=0,4)
        enddo
        if(ubi=='psi'.or.ubi=='bsp') then
          do i_omega=0,n_omega
            write(4,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(sv(i_omega,i))*evolt/2,i=0,4)
          enddo
        endif
        close(3)
        if(ubi=='psi'.or.ubi=='bsp') close(4)
c ------------ At Kf -----------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_at_Kf.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_at_Kf.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_at_Kf.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_at_Kf.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=n_k_fermi-2,n_k_fermi+2)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_at_Kf.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_at_Kf.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_at_Kf.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_at_Kf.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=n_k_fermi-2,n_k_fermi+2)
        enddo
        close(3)
c ------------ Large K -----------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_large_K.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_large_K.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_large_K.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Re_large_K.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=n_k_heg-4,n_k_heg)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_large_K.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_large_K.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_large_K.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIGC_omega_Im_large_K.bsp')
        endif
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=n_k_heg-4,n_k_heg)
        enddo
        close(3)
c --------------------------------------------------------------        
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Re.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Re.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Re.psi')
          open(4,file=allfile(1:i_len)//'_SIG_VRT_K_Re.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Re.bsp')
          open(4,file=allfile(1:i_len)//'_SIG_VRT_K_Re.bsp')
        endif
        do i=0,n_k_heg
          write(3,'(f9.3,5(1x,f12.7))')k_m_heg(i),
     &   ((sigx_k_heg(i,nspin)+real(s(i_omega,i)))*evolt/2,i_omega=0,4)
        enddo
        if(ubi=='psi'.or.ubi=='bsp') then
          do i=0,n_k_heg
            write(4,'(f9.3,5(1x,f12.7))')k_m_heg(i),
     &         (real(sv(i_omega,i))*evolt/2,i_omega=0,4)
          enddo
        endif
        close(3)
        if(ubi=='psi'.or.ubi=='bsp') close(4)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Im.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Im.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Im.psi')
          open(4,file=allfile(1:i_len)//'_SIG_VRT_K_Im.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_k_Im.bsp')
          open(4,file=allfile(1:i_len)//'_SIG_VRT_K_Im.bsp')
        endif
        do i=0,n_k_heg
          write(3,'(f9.3,5(1x,f12.7))')k_m_heg(i),
     &    (imag(s(i_omega,i))*evolt/2,i_omega=0,4)
        enddo
        if(ubi=='psi'.or.ubi=='bsp') then
          do i=0,n_k_heg
            write(4,'(f9.3,5(1x,f12.7))')k_m_heg(i),
     &               (imag(sv(i_omega,i))*evolt/2,i_omega=0,4)
          enddo
        endif
        close(3)
        if(ubi=='psi'.or.ubi=='bsp') close(4)
      endif
      deallocate(s,sv)
      end
      
      
      
      
      subroutine output_heg1_p_nu(pnu,key)
c     key = 0 for P
c     key = 1 for P_WI
c     key = 2 for P_0
c     key = 3 for P_corrected       
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
      integer, intent(in) :: key
      real*8, intent(in) :: pnu(0:n_q_heg,0:n_nu)
      character*6 ch
      integer :: i_len,i_nu,i,n,n1,n2
      real*8, allocatable :: s(:,:)
      ch(1:2)='_P'
      n=2
      if(key==1) then
        ch='_P_WI'
        n=5
      endif
      if(key==2) then
        ch='_P0'
        n=3
      endif
      if(key==3) then
        ch='_Pcorr'
        n=6
      endif
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=pnu(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_nu.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_nu.dft')
        endif
        write(3,*)'Output at Q: 0,1,n_q/4,n_q/2,n_q'
        n1=n_q_heg/4
        n2=n_q_heg/2
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,1),s(i_nu,n1),s(i_nu,n2),
     &                                 s(i_nu,n_q_heg)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q.dft')
        endif
        write(3,*)'Output at NU: 0,1,2,n_exa,n_nu'
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                             (s(i_nu,i),i_nu=0,2),
     &     s(n_nu_exa,i),s(n_nu,i)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg_p_tau(ptau,key)
c     key = 0 for P
c     key = 1 for P_0      
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
      integer, intent(in) :: key
      real*8, intent(in) :: ptau(0:n_q_heg,0:n_tau/2)
      character*5 ch
      integer :: i_len,i_tau,i,n,n1,n2
      real*8, allocatable :: s(:,:)
      ch(1:2)='_P'
      n=2
      if(key==1) then
        ch='_P0'
        n=3
      endif
      allocate(s(0:n_tau/2,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=ptau(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_tau.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_tau.dft')
        endif
        write(3,*)'Output at Q: 0,1,n_q/4,n_q/2,n_q'
        n1=n_q_heg/4
        n2=n_q_heg/2
        do i_tau=0,n_tau/2
          write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &	                   (s(i_tau,i),i=0,1),s(i_tau,n1),s(i_tau,n2),
     &                                 s(i_tau,n_q_heg)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q_tau.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//ch(1:n)//'_Q_tau.dft')
        endif
        write(3,*)'Output at TAU: 0,1,2,n_tau/2-1,n_tau/2'
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                             (s(i_tau,i),i_tau=0,2),
     &     s(n_tau/2-1,i),s(n_tau/2,i)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      subroutine output_heg1_w_nu
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
      integer :: i_len,i_nu,i,ir
      real*8 :: r1
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=w_q_nu_heg(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- W (nu;q) -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_W_nu_q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_W_nu_q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_W_nu_q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_W_nu_q.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                  (s(i_nu,i),i=0,4)
        enddo
        close(3)
c -------- W (q;nu) -------------------------
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_W_Q_Nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_W_Q_Nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_W_Q_Nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_W_Q_Nu.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,e12.5))')q_m_heg(i),v_q_heg(i),
     &	                               (s(i_nu,i),i_nu=0,3)
        enddo
        close(3)
c -------- W (r;nu) -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_W_r_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_W_r_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_W_r_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_W_r_nu.bsp')
        endif
        do ir=0,n_r_heg
          r1=r_m_heg(ir)
          write(3,'(f9.3,6(1x,f9.4))')r1,
     &	       (w_r_nu_heg(ir,i_nu),i_nu=0,5)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      subroutine output_heg1_g_rax
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
      integer :: i_len,i,i_omega,i0,ik,ispin
      real*8 :: deltae,w,a(0:4),eqp,d,d0,f0,fm,fp,aaa,bbb,ccc
      real*8, allocatable :: aa(:)
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_A.gw')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_A.dft')
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_A.hf')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_A.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_A.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_A.bsp')
        endif
      endif
      if(maswrk) then
        write(3,*)
     & ' First column - frequency(eV), second - Full Spectral function'
        write(3,*)
     &	' Further columns - Spectral functions for first few K'
        deltae=(emaxdos-emindos)/ndos
        do i_omega=0,ndos
          w=(emindos+i_omega*deltae)*evolt/2.d0
          a(0)=0.d0
          do i=0,n_k_heg
c            a(0)=a(0)-sum(imag(g_rax_heg(i_omega,i,:)))/pi*wgt(i)
          enddo
          do i=1,4
            i0=i-1
            a(i)=-sum(imag(g_rax_heg(i_omega,i0,:)))/pi
          enddo
          if(nspin==1) a=a*2.d0    !! spin
          write(3,'(f9.3,5(1x,f12.7))')w,(a(i),i=0,4)
        enddo
      endif
      if(maswrk) close(3)
c ---------- Spectrum -------------------------------------
      spectrum_ac=1.d5
      allocate(aa(0:ndos))
      do ik=0,n_k_heg
        aa=0.d0
        do ispin=1,nspin
          eqp=k_m_heg(ik)**2+vxc_heg(ispin)-chem_pot_lda
          d=1.d5
          i0=-1
          do i_omega=0,ndos
            w=emindos+i_omega*deltae
            aa(i_omega)=aa(i_omega)-imag(g_rax_heg(i_omega,ik,ispin))
     &                             /pi
            d0=abs(eqp-w)
            if(d0<d) then
              d=d0
              i0=i_omega
            endif
          enddo
        enddo
        if(nspin==1) aa=aa*2.d0    !! spin
        if(i0<=0.or.i0>ndos-1) cycle
        if(aa(i0-1)<aa(i0).and.aa(i0)<aa(i0+1)) then
          do i=i0+1,ndos-1
            if(aa(i-1)<aa(i).and.aa(i)>aa(i+1)) then
              fm=aa(i-1)
              f0=aa(i)
              fp=aa(i+1)
              ccc=f0
              bbb=(fp-fm)/2/deltae
              aaa=(fp-bbb*deltae-ccc)/deltae**2
              spectrum_ac(ik)=emindos+i*deltae-bbb/2/aaa
              exit
            endif
          enddo
        else if(aa(i0-1)<aa(i0).and.aa(i0)>aa(i0+1)) then
          fm=aa(i0-1)
          f0=aa(i0)
          fp=aa(i0+1)
          ccc=f0
          bbb=(fp-fm)/2/deltae
          aaa=(fp-bbb*deltae-ccc)/deltae**2
          spectrum_ac(ik)=emindos+i0*deltae-bbb/2/aaa
        else if(aa(i0-1)>aa(i0).and.aa(i0)>aa(i0+1)) then
          do i=i0-1,1,-1
            if(aa(i-1)<aa(i).and.aa(i)>aa(i+1)) then
              fm=aa(i-1)
              f0=aa(i)
              fp=aa(i+1)
              ccc=f0
              bbb=(fp-fm)/2/deltae
              aaa=(fp-bbb*deltae-ccc)/deltae**2
              spectrum_ac(ik)=emindos+i*deltae-bbb/2/aaa
              exit
            endif
          enddo
        endif
      enddo
      deallocate(aa)
      end
      
      
      
      
      subroutine output_heg1_g_pair  
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
      integer :: i_len,ir,i
      real*8 :: r1,vr
      if(maswrk) then
        i_len=len_trim(allfile)  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_g_pair.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_g_pair.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_g_pair.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_g_pair.bsp')
        endif
        do ir=0,n_r_heg
          r1=r_m_heg(ir)
          vr=1.d0-resp_r_tau_heg(ir,0)/ro_el_gas**2
          write(3,'(f12.6,5(1x,f12.7))')r1,vr
        enddo
        close(3) 
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SSF.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SSF.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SSF.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SSF.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f12.6,5(1x,f12.7))')q_m_heg(i),
     &	                            -resp_q_tau_heg(i,0)/ro_el_gas
        enddo
        close(3)
      endif
      end
      
      
      
      
      subroutine output_heg1_e_nu
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
      integer :: i_len,i_nu,i
      real*8, allocatable :: y(:,:),s(:,:)
      allocate(s(0:n_nu,0:n_q_heg))
      allocate(y(0:n_q_heg,0:n_nu))
      do i=0,n_q_heg
        s(:,i)=edd_q_nu_heg(i,:)
        y(i,:)=e1_heg(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Density-density Dielectric function ----------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_E_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,4)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_E_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E_Q.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                               (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
c ---- Eps^(-1) ------------------------
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_E1_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E1_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E1_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E1_Q.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                               (y(i,i_nu),i_nu=0,4)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_E1_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E1_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E1_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E1_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (y(i,i_nu),i=0,4)
        enddo
        close(3)  
      endif
      deallocate(s,y)
      end
      
      
      
      
      
      
      subroutine output_heg1_wi(wi_chk,n_chk)
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
      integer, intent(in) :: n_chk(4)
      real*8, intent(in) :: wi_chk(4,2)
      integer :: i
      if(.not.allocated(num_omega1_adapt)) return
      if(nproc_k/=1) then
        call dgop(wi_chk(1,2),4,'  +',comm_k)
        call dgop(wi_chk(1,1),4,'max',comm_k)
        call igop(n_chk,1,'  +',comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(wi_chk(1,2),4,'  +',comm_t)
        call dgop(wi_chk(1,1),4,'max',comm_t)
        call igop(n_chk,1,'  +',comm_t)
      endif
      if(maswrk) then
        write(iun,*)' WI checking: q0nu0 nu0  q0 full:'
        write(iun,'(a3,a9,4(1x,e12.5))')ubi,' Average:',
     &	                                (wi_chk(i,2)/n_chk(i),i=1,4)
        write(iun,'(a3,a9,4(1x,e12.5))')ubi,' Maximal:',
     &	                                (wi_chk(i,1),i=1,4)
      endif
      end
      
      
      subroutine output_heg1_z_ren
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
      integer :: i_len,i
      character*3 :: ch
      real*8 :: a
      if(ubi==' qp') ch='qp '
      if(ubi==' gw') ch='gw '
      if(ubi=='psi') ch='psi'
      if(ubi=='bsp') ch='bsp'
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_Z_REN_K.'//ch)
        do i=0,n_k_heg
          write(3,'(f9.4,1x,f12.7)')k_m_heg(i),1.d0
     &	                            /z_ren_heg(i,nspin)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_Z_0_K.'//ch)
        do i=0,n_k_heg
          a=sig_0_heg(i,nspin)+sigx_k_heg(i,nspin)
          write(3,'(f9.4,1x,f12.7)')k_m_heg(i),a
        enddo
        close(3)
      endif
      end
      
      
      
      
      subroutine output_sig_heg1_real
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
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: i_len,i,i_omega,j,k0,k
      real*8 :: deltae,w,w0,w2,w4,e0,e1,e2,d,d1,a,b,c,pi2a,eqp
      real*8, allocatable :: e(:,:)
      pi2a=(pi+pi)/par
      allocate(e(0:ndos,0:n_k_heg))
      spectrum=0.d0
      deltae=(emaxdos-emindos)/ndos
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi/='dft'.and.ubi/=' hf') then
          if(ubi==' gw') then
            open(3,file=allfile(1:i_len)//'_Sig_real.gw')
          else if(ubi==' qp') then
            open(3,file=allfile(1:i_len)//'_Sig_real.qp')
          else if(ubi=='psi') then
            open(3,file=allfile(1:i_len)//'_Sig_real.psi')
          else if(ubi=='bsp') then
            open(3,file=allfile(1:i_len)//'_Sig_real.bsp')
          endif
          do i_omega=0,ndos
            w=emindos+i_omega*deltae
            do i=0,n_k_heg
              e(i_omega,i)=e_star_x(i,nspin)
     &                    +real(sig_real_heg(i,i_omega,nspin))
     &                    -chem_pot-w
            enddo
            w=w*evolt/2.d0
            write(3,'(f9.3,5(1x,f12.7))')w,
     &           ((real(sig_real_heg(i,i_omega,nspin))
     &             +sigx_k_heg(i,nspin))*evolt/2,
     &                             i=0,4)
          enddo
          close(3)
          if(ubi==' gw') then
            open(3,file=allfile(1:i_len)//'_Spectrum_Omega.gw')
          else if(ubi==' qp') then
            open(3,file=allfile(1:i_len)//'_Spectrum_Omega.qp')
          else if(ubi=='psi') then
            open(3,file=allfile(1:i_len)//'_Spectrum_Omega.psi')
          else if(ubi=='bsp') then
            open(3,file=allfile(1:i_len)//'_Spectrum_Omega.bsp')
          endif
          do i_omega=0,ndos
            w=(emindos+i_omega*deltae)*evolt/2.d0
            write(3,'(f9.3,5(1x,f12.7))')w,(e(i_omega,i)*evolt/2,
     &	                                    i=0,4)
          enddo
          close(3)
        endif
c ------ One electron spectrum ------------------------------
        w0=5
        w2=10.d0*deltae**2
        w4=34.d0*deltae**4
        spectrum=10000.d0
        do i=0,n_k_heg
c ------ First we locate the quasiparticle energy ------------
          eqp=e_star_qp(i,nspin)-chem_pot
          if(ubi==' qp') then
            spectrum(i)=eqp
          else if(ubi=='dft') then
            spectrum(i)=e_star_free(i,nspin)+vxc_heg(nspin)-chem_pot
          else if(ubi==' hf') then
            spectrum(i)=e_star_free(i,nspin)+sigx_k_heg(i,nspin)
     &                                       -chem_pot
          else
c ------ Now we locate the true energy closest to E_qp -----------
            j=0
            d=abs(emindos-eqp)
            do i_omega=0,ndos-1
              if(e(i_omega,i)*e(i_omega+1,i)<=0.d0) then
                w=emindos+i_omega*deltae
                d1=abs(w-eqp)
                if(d1<d) then
                  j=i_omega
                  d=d1
                endif
              endif
            enddo
            if(j<2.or.j>ndos-2) cycle
            e0=0.d0
            e1=0.d0
            e2=0.d0
            do k0=-2,2
              k=j+k0
              w=k0*deltae
              e0=e0+e(k,i)
              e1=e1+w*e(k,i)
              e2=e2+w*w*e(k,i)
            enddo
            b=e1/w2
            c=(e0-w0*e2/w2)/(w2-w0*w4/w2)
            a=(e0-w2*c)/w0
            d=sqrt(b*b-4.d0*a*c)
            w=(-b-d)/2/c
            if(w>-deltae.and.w<deltae) then
              spectrum(i)=w+deltae*j+emindos
            else
              spectrum(i)=(-b+d)/2/c+deltae*j+emindos
            endif
          endif
        enddo   !! over i
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.dft')
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_Spectrum_K.hf')
        endif
        if(ubi/='dft'.and.ubi/=' hf') then
          write(3,*)' Second column - spectrum from zeroes of denom-r'
          write(3,*)' Third column - spectrum from peaks of spectr. f-n'
        endif
        do i=0,n_k_heg
          if(ubi/='dft'.and.ubi/=' hf') then
            write(3,'(f9.4,2(1x,f12.7))')k_m_heg(i),
     & 	                                 spectrum(i)*evolt/2,
     & 	                                 spectrum_ac(i)*evolt/2
          else
            write(3,'(f9.4,1x,f12.7)')k_m_heg(i),
     & 	                                 spectrum(i)*evolt/2
          endif
        enddo
        close(3)            
      endif
      deallocate(e)
#ifdef MPI
      if(goparr) call brdcst(msgdbl,spectrum,8*(n_k_heg+1),master,
     & 	                   MPI_COMM_WORLD)
#endif     
      end
      
      
      
      subroutine output_heg1_e1_rax
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: i,nrax
      real*8 :: h_xi_rax
      real*8, allocatable :: wrx(:)
      complex*16, allocatable :: e1(:)
      allocate(e1(0:n_nu))
      nrax=nrax_chi
      allocate(wrx(0:nrax))
      h_xi_rax=sqrt(freq_chi)/nrax
      do i=0,nrax
        wrx(i)=(i*h_xi_rax)**2
      enddo
      do i=1,nstar
        e1=e1_heg(:,i)-1.d0
        call pade_boson_all(e1,wrx,nrax,e1_rax_heg(0,1,i),2)
        call pade_boson_all(e1,-wrx,nrax,e1_rax_heg(0,2,i),2)
      enddo  !! over i
      e1_rax_heg=e1_rax_heg+1.d0
      deallocate(wrx,e1)
      end
      
      
      
      
      subroutine output_heg1_resp_nu(key)
c     key = 0 - Noninteracting (DFT) response
c     key = 1 - Interacting (enhanced) response       
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
      integer, intent(in) :: key
      character*1 :: ch
      integer :: i_len,i_nu,i
      real*8, allocatable :: s(:,:)
      if(key==0) ch='0'
      if(key==1) ch='_'
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        if(key==0) s(:,i)=resp0_q_nu_heg(i,:)
        if(key==1) s(:,i)=resp_q_nu_heg(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_nu.dft')	  
        else if(ubi==' gw') then	 
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,4)
        enddo
        close(3)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_Q.dft')	  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP'//ch//'_Q.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg1_f_xc
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
      integer :: i_len,i_nu,i
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=f_xc(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------  
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_F_XC_nu.dft')	  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_F_XC_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_F_XC_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_F_XC_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_F_XC_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,4)
        enddo
        close(3) 
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_F_XC_Q.dft')	  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_F_XC_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_F_XC_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_F_XC_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_F_XC_Q.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg1_p_nu_vrt
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
      integer :: i_len,i_nu,i
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=p_q_nu_vertex_rad(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c ------- Polarizability -------------------------	  
        if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_P_nu_VRT.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_P_nu_VRT.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_P_nu_VRT.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,4)
        enddo
        close(3)
        if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_P_Q_VRT.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_P_Q_VRT.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_P_Q_VRT.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                               (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end 
      
      
      
      
      
      subroutine output_heg1_g_fac
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
      integer :: i_len,i_nu,i
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=g_fac_s(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------  
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_nu.dft')	  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (s(i_nu,i),i=0,4)
        enddo
        close(3) 
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_Q.dft')	  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_G_FAC_S_Q.bsp')
        endif
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      subroutine output_heg1_sig_x
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
      integer :: i_len,i
      real*8 :: pi2a,k,y,sm,sa1,pi4
      pi4=4.d0*pi
      pi2a=(pi+pi)/par
      if(maswrk) then
        i_len=len_trim(allfile)
c ----------- Exchange ------------------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_SIG_X_k_Re.gw')
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_SIG_X_k_Re.hf')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_SIG_X_k_Re.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_SIG_X_k_Re.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_SIG_X_k_Re.bsp')
        endif
        write(3,'(a37)')'   K    Sig_X(K)   Sig_x_asymp  Mahan'
        do i=0,n_k_heg
          k=k_m_heg(i)
          y=k/k_fermi
          if(i==0) then
            sm=-2.d0
          else if(abs(y-1.d0)<1.d-12) then
            sm=-1.d0
          else
            sm=abs((1.d0+y)/(1.d0-y))
            sm=(1.d0-y*y)*log(sm)/2.d0/y
            sm=-(1.d0+sm)
          endif
          sm=2.d0*k_fermi*sm/pi
          sa1=0.d0
          if(i/=0) sa1=-pi4*ro_el_gas/k**2
          write(3,'(f9.3,3(1x,f12.7))')k,sigx_k_heg(i,nspin)*evolt/2,
     &	                             sa1*evolt/2,sm*evolt/2
        enddo
        close(3)
      endif
      end
      
      
      
      
      subroutine output_heg1_resp_line
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
      integer :: i_len,i_nu,i,i0
      real*8 :: pi2a,v
      real*8, allocatable :: e1(:,:)
      pi2a=(pi+pi)/par
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Charge Response -------------------------
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.dft')	  
        else if(ubi==' hf') then	 
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.hf') 
        else if(ubi==' gw') then	 
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (resp_line_heg(i,i_nu),i=1,min(5,n_line))
        enddo
        close(3)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.dft')	  
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.hf')  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Q.bsp')
        endif
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (resp_line_heg(i,i_nu),i_nu=0,4)
        enddo
        close(3)
c -------- Spin ZZ Response -------------------------
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.dft')	  
        else if(ubi==' hf') then	 
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.hf') 
        else if(ubi==' gw') then	 
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (resp_zz_line_heg(i,i_nu),i=1,min(5,n_line))
        enddo
        close(3)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.dft')	  
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.hf')  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP_ZZ_!!!_Q.bsp')
        endif
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (resp_zz_line_heg(i,i_nu),i_nu=0,4)
        enddo
        close(3)
c -------- Inverse dielectric function -------------------------
        allocate(e1(n_line,0:n_nu))
        do i=1,n_line
          i0=kline_in_npnt(i)
          v=v_q_heg(i0)
          do i_nu=0,n_nu
            e1(i,i_nu)=1.d0+v*resp_line_heg(i,i_nu)
          enddo
        enddo
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.dft')	  
        else if(ubi==' hf') then	 
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.hf') 
        else if(ubi==' gw') then	 
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_nu.bsp')
        endif
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (e1(i,i_nu),i=1,min(5,n_line))
        enddo
        close(3)
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.dft')	  
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.hf')  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_E1_!!!_Q.bsp')
        endif
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (e1(i,i_nu),i_nu=0,4)
        enddo
        close(3)
        deallocate(e1)
      endif
      end
      
      
      
      
      subroutine output_heg1_resp_static
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
      integer :: i_len,i,i0
      real*8 :: pi2a,v
      real*8, allocatable :: e1(:)
      pi2a=(pi+pi)/par
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Inverse dielectric function -------------------------
        allocate(e1(n_line))
        do i=1,n_line
          i0=kline_in_npnt(i)
          v=v_q_heg(i0)
          e1(i)=1.d0+v*resp_line_heg(i,0)
        enddo
c -------- Static Responses -------------------------
        if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.dft')	  
        else if(ubi==' hf') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.hf')  
        else if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_RESP_!!!_Static.bsp')
        endif
        write(3,*)' Q       R0      R_Charge    R_ZZ     E1 '
        do i=1,n_line
          i0=kline_in_npnt(i)
          write(3,'(f9.3,4(1x,f12.7))')pi2a*qabs_line(i),
     &	                   resp0_q_nu_heg(0,i0),resp_line_heg(i,0),
     &                       resp_zz_line_heg(i,0),e1(i)
        enddo
        close(3)
        deallocate(e1)
      endif
      end
      
      
      
      
      subroutine output_heg1_v_r(key)
c     key = 0  - only V(r)_exact | V(r)
c     key = 1  - V(r)_exact | V(r)  | W_scr(r,tau=0) | W_full(r,tau=0) 
c                                   | W_scr(r,nu=0)  | W_full(r,nu=0)     
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
      integer, intent(in) :: key
      integer :: i_len,ir
      real*8 :: r1,vr
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_V_r')
        do ir=0,n_r_heg
          r1=r_m_heg(ir)
          vr=0.d0
          if(ir/=0) vr=2.d0/r1
          if(key==0) write(3,'(f9.3,5(1x,f12.7))')r1,vr,v_r_heg(ir)
          if(key==1) write(3,'(f9.3,6(1x,f9.4))')r1,vr,v_r_heg(ir),
     &	    w_r_tau_heg(ir,0),v_r_heg(ir)+w_r_tau_heg(ir,0),
     &	    w_r_nu_heg(ir,0),v_r_heg(ir)+w_r_nu_heg(ir,0)
        enddo
        close(3)
      endif
      end
      
       
      
      
      
      subroutine output_heg_sig0_omega
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
      integer :: i_len,i_omega,i
      complex*16, allocatable :: s(:,:)
      allocate(s(0:n_omega,0:n_k_heg))
      do i=0,n_k_heg
        s(:,i)=sigc0_k_omega_heg(i,:,nspin)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c ------------ Small K ----------------------------------------- 
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Re.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Im.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=0,4)
        enddo
        close(3)
c ------------ At Kf -----------------------------------------
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Re_at_Kf.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=n_k_fermi-2,n_k_fermi+2)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Im_at_Kf.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=n_k_fermi-2,n_k_fermi+2)
        enddo
        close(3)
c ------------ Large K -----------------------------------------
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Re_large_K.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &    (real(s(i_omega,i))*evolt/2,i=n_k_heg-4,n_k_heg)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_SIGC0_omega_Im_large_K.gw')
        do i_omega=0,n_omega
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &            (imag(s(i_omega,i))*evolt/2,i=n_k_heg-4,n_k_heg)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg_w_tau
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
      integer :: i_len,i_tau,i,n1,n2
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_tau/2,0:n_q_heg))
      do i=0,n_q_heg
        s(:,i)=w_q_tau_heg(i,:)
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_W_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_W_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_W_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_W_tau.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_W_tau.dft')
        endif
        write(3,*)'Output at Q: 0,1,n_q/4,n_q/2,n_q'
        n1=n_q_heg/4
        n2=n_q_heg/2
        do i_tau=0,n_tau/2
          write(3,'(f9.3,2(1x,f16.4),3(1x,f12.7))')tau_mesh(i_tau),
     &	                   (s(i_tau,i),i=0,1),s(i_tau,n1),s(i_tau,n2),
     &                                 s(i_tau,n_q_heg)
        enddo
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_W_Q_tau.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_W_Q_tau.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_W_Q_tau.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_W_Q_tau.bsp')
        else if(ubi=='dft') then
          open(3,file=allfile(1:i_len)//'_W_Q_tau.dft')
        endif
        write(3,*)'Output at TAU: 0,1,2,n_tau/2-1,n_tau/2'
        do i=0,n_q_heg
          write(3,'(f9.3,5(1x,f12.7))')q_m_heg(i),
     &	                             (s(i_tau,i),i_tau=0,2),
     &     s(n_tau/2-1,i),s(n_tau/2,i)
        enddo
        close(3)
      endif
      deallocate(s)
      end