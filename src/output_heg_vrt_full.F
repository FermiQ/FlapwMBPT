      subroutine output_heg_vrt_full(vrt,vrt_k_q,vrt_nu_q,vrt_om_nu,
     &                               vrt_om_q,vrt_om_k,vrt_nu_k)
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
      complex*16, intent(in) :: vrt(nc_line,0:n_nu),
     &                          vrt_k_q(nc_line,nc_line),
     &                          vrt_nu_q(0:n_nu_exa,nc_line),
     &                          vrt_om_nu(0:n_omega_exa,0:n_nu_exa),
     &                          vrt_om_q(0:n_omega_exa,nc_line),
     &                          vrt_om_k(0:n_omega_exa,nc_line),
     &                          vrt_nu_k(0:n_nu_exa,nc_line)
      integer :: i_len,i_nu,i,iq,ik,i_omega
      real*8 :: pi2a
      pi2a=(pi+pi)/par
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_nu_Re.'//ubi)
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (real(vrt(i,i_nu)),i=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_nu_Im.'//ubi)
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (imag(vrt(i,i_nu)),i=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_Q_Re.'//ubi)
        do i=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(i),
     &	                   (real(vrt(i,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_Q_Im.'//ubi)
        do i=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(i),
     &	                   (imag(vrt(i,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_K_Q_Re.'//ubi)
        do ik=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(ik),
     &	                   (real(vrt_k_q(ik,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_K_Q_Im.'//ubi)
        do ik=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(ik),
     &	                   (imag(vrt_k_q(ik,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_Q_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_nu_q(i_nu,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_Q_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_nu_q(i_nu,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_NU_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_nu(i_omega,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_NU_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_nu(i_omega,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_OM_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_om_nu(i_omega,i_nu)),i_omega=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_OM_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_om_nu(i_omega,i_nu)),i_omega=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_Q_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_q(i_omega,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_Q_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_q(i_omega,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_K_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_k(i_omega,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_K_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_k(i_omega,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_K_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_nu_k(i_nu,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_K_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_nu_k(i_nu,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
      endif
      end
      
      
      
      
      subroutine output_heg_vrt_spl(vrt,vrt_k_q,vrt_nu_q,vrt_om_nu,
     &                               vrt_om_q,vrt_om_k,vrt_nu_k)
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
      complex*16, intent(in) :: vrt(nc_line,0:n_nu),
     &                          vrt_k_q(nc_line,nc_line),
     &                          vrt_nu_q(0:n_nu_exa,nc_line),
     &                          vrt_om_nu(0:n_omega_exa,0:n_nu_exa),
     &                          vrt_om_q(0:n_omega_exa,nc_line),
     &                          vrt_om_k(0:n_omega_exa,nc_line),
     &                          vrt_nu_k(0:n_nu_exa,nc_line)
      integer :: i_len,i_nu,i,iq,ik,i_omega
      real*8 :: pi2a,v1,v2
      pi2a=(pi+pi)/par
      if(maswrk) then
        write(iun,*)' Q-limit of the vertex is ',real(vrt(1,0))
        v1=real(vrt(1,1))
        v2=real(vrt(1,2))
        write(iun,*)' Nu-limit of the vertex is ',(4*v1-v2)/3
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_nu_Re.'//ubi)
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (real(vrt(i,i_nu)),i=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_nu_Im.'//ubi)
        do i_nu=0,n_nu
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &	                   (imag(vrt(i,i_nu)),i=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_Q_Re.'//ubi)
        do i=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(i),
     &	                   (real(vrt(i,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_FULL_Q_Im.'//ubi)
        do i=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(i),
     &	                   (imag(vrt(i,i_nu)),i_nu=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_K_Q_Re.'//ubi)
        do ik=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(ik),
     &	                   (real(vrt_k_q(ik,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_K_Q_Im.'//ubi)
        do ik=1,nc_line
          write(3,'(f9.3,5(1x,f12.7))')pi2a*qabs_line_c(ik),
     &	                   (imag(vrt_k_q(ik,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_Q_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_nu_q(i_nu,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_Q_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_nu_q(i_nu,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_NU_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_nu(i_omega,i_nu)),
     &                        i_nu=0,min(4,n_nu_exa))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_NU_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_nu(i_omega,i_nu)),
     &                        i_nu=0,min(4,n_nu_exa))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_OM_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_om_nu(i_omega,i_nu)),
     &                  i_omega=0,min(4,n_omega_exa))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_OM_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_om_nu(i_omega,i_nu)),
     &              i_omega=0,min(4,n_nu_exa))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_Q_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_q(i_omega,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_Q_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_q(i_omega,iq)),iq=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_K_Re.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (real(vrt_om_k(i_omega,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_OM_K_Im.'//ubi)
        do i_omega=0,n_omega_exa
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                (imag(vrt_om_k(i_omega,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_K_Re.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (real(vrt_nu_k(i_nu,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_NU_K_Im.'//ubi)
        do i_nu=0,n_nu_exa
          write(3,'(f9.3,5(1x,f12.7))')w_nu(i_nu)*evolt/2,
     &                (imag(vrt_nu_k(i_nu,ik)),ik=1,min(5,nc_line))
        enddo
        close(3)
      endif
      end