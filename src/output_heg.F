      subroutine output_heg_resp_nu(key)
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
      integer :: i_len,i_nu,i,i0
      real*8 :: pi2a
      real*8, allocatable :: s(:,:)
      if(key==0) ch='0'
      if(key==1) ch='_'
      allocate(s(0:n_nu,n_line))
      do i=1,n_line
        i0=kline_in_npnt(i)
        if(key==0) s(:,i)=resp0_q_nu_heg(:,i0)
        if(key==1) s(:,i)=resp_q_nu_heg(:,i0)
      enddo
      pi2a=(pi+pi)/par
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
     &	                   (s(i_nu,i),i=1,min(5,n_line))
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
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg_f_xc
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
      real*8 :: pi2a
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,n_line))
      do i=1,n_line
        i0=kline_in_npnt(i)
        s(:,i)=f_xc(i0,:)
      enddo
      pi2a=(pi+pi)/par
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
     &	                   (s(i_nu,i),i=1,min(5,n_line))
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
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_ladder_conv(key)
c     key = 0 - For Nu_c - mesh
c     key = 1 - For Nu_d_c - mesh            
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
      integer :: i_len,i_nu,i
      real*8 :: pi2a
      real*8, allocatable :: s(:,:)
      if(key==0) allocate(s(0:n_nu,nc_line))
      if(key==1) allocate(s(0:n_nu,nc_line))
      do i=1,nc_line
        s(:,i)=conv_ladder(i,:)
      enddo
      pi2a=(pi+pi)/par
      if(maswrk) then
        i_len=len_trim(allfile)
c -------- Polarizability -------------------------	  
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_nu.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_nu.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_nu.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_nu.bsp')
        endif
        if(key==0) then
          do i_nu=0,n_nu
            write(3,'(f9.3,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                     (s(i_nu,i),i=1,min(5,nc_line))
          enddo
        else if(key==1) then
          do i_nu=0,n_nu
            write(3,'(f9.3,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                     (s(i_nu,i),i=1,min(5,nc_line))
          enddo
        endif
        close(3)
        if(ubi==' gw') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_Q.gw')
        else if(ubi==' qp') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_Q.qp')
        else if(ubi=='psi') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_Q.psi')
        else if(ubi=='bsp') then
          open(3,file=allfile(1:i_len)//'_Conv_ladder_Q.bsp')
        endif
        do i=1,nc_line
          write(3,'(f9.3,5(1x,e12.5))')pi2a*qabs_line_c(i),
     &	                               (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      
      subroutine output_heg_g_fac
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
      real*8 :: pi2a
      real*8, allocatable :: s(:,:)
      allocate(s(0:n_nu,n_line))
      do i=1,n_line
        i0=kline_in_npnt(i)
        s(:,i)=g_fac_s(i0,:)
      enddo
      pi2a=(pi+pi)/par
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
     &	                   (s(i_nu,i),i=1,min(5,n_line))
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
        do i=1,n_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line(i),
     &	                   (s(i_nu,i),i_nu=0,4)
        enddo
        close(3)
      endif
      deallocate(s)
      end
      
      
      
      subroutine output_heg_sig_x
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
        do i=1,n_line
          i0=kline_in_npnt(i)
          k=pi2a*qabs_star(i0)
          y=k/k_fermi
          if(i0==1) then
            sm=-2.d0
          else if(abs(y-1.d0)<1.d-12) then
            sm=-1.d0
          else
            sm=abs((1.d0+y)/(1.d0-y))
            sm=(1.d0-y*y)*log(sm)/2.d0/y
            sm=-(1.d0+sm)
          endif
          sm=2.d0*k_fermi*sm/pi
          sa1=-pi4*ro_el_gas/k**2
          write(3,'(f9.3,3(1x,f12.7))')k,sigx_k_heg(i0,nspin)*evolt/2,
     &	                             sa1*evolt/2,sm*evolt/2
        enddo
        close(3)
      endif
      end
      
      
      
      
      subroutine output_heg_resp_line
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
      
      
      
      
      subroutine output_heg_resp_static
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