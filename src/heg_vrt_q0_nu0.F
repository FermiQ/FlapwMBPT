	subroutine heg_vrt_q0_nu0(wqt)
	use heg_mod
	use manager_mod
      use parallel_mod
	use solid_mod
	use units_mod
	implicit none
      real*8, intent(in) :: wqt(0:n_q_heg,0:n_tau/2)
	integer :: ik,ispin,iq,ix,i,i_tau,i0_tau,i_omega,i_len
      real*8 :: k,q,k2q2,kq2,x,p,k2,p2,q2,hh,pmi,pma,ymi,yma,y,ax0,
     &          aa(0:4)
      real*8, allocatable :: g(:),ax(:),xm(:),xw(:),
     &                       vrt_q0_nu0_heg(:,:,:),g2_tau(:,:),
     &                       g2_spl(:,:,:),
     &                       vrt_stat_q0_nu0_heg(:,:)
      complex*16, allocatable :: g2_om(:),vrt_om_q0_nu0_heg(:,:,:)
      allocate(vrt_om_q0_nu0_heg(0:n_k_heg,0:n_omega,nspin))
      allocate(vrt_stat_q0_nu0_heg(0:n_k_heg,nspin))
      allocate(vrt_q0_nu0_heg(0:n_k_heg,0:n_tau,nspin))
      allocate(g(0:n_tau))
      allocate(ax(0:n_tau))
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      allocate(g2_om(0:n_omega))
      allocate(g2_tau(0:n_k_heg,0:n_tau))
      allocate(g2_spl(n_k_heg,4,0:n_tau))
      vrt_q0_nu0_heg=0.d0
      vrt_stat_q0_nu0_heg=0.d0
	do ispin=1,nspin
	  do ik=0,n_k_heg
          do i_omega=0,n_omega
            g2_om(i_omega)=-g_k_omega_heg(ik,i_omega,ispin)**2
          enddo
c          call from_omega_to_tau_v(1,g2_om,g2_tau(ik,:))
        enddo
c -------- Splines for G2(k,tau) --------------------------------  
        do i_tau=0,n_tau
          call spline_inhmg(k_m_heg,g2_tau(0,i_tau),g2_spl(1,1,i_tau),
     &                      g2_spl(1,2,i_tau),g2_spl(1,3,i_tau),
     &                      g2_spl(1,4,i_tau),n_k_heg,0,0.d0,0.d0)
        enddo
c --------------------------------------------------------------        
	  do ik=0,n_k_heg
          k=k_m_heg(ik)
          k2=k*k
          do iq=0,n_q_heg
            q=q_m_heg(iq)
            q2=q*q
            k2q2=k2+q2
            kq2=2.d0*k*q
c -------- Adaptive X-mesh ----------------------------------
            pmi=abs(k-q)
            pma=k+q
            ymi=x_fermi+asinh((pmi-k_fermi)/alf_k_heg)
            yma=x_fermi+asinh((pma-k_fermi)/alf_k_heg)
            hh=(yma-ymi)/n_x_heg
            xm=0.d0
            xw=2.d0/(n_x_heg+1)
            if(k>1.d-8.and.q>1.d-8) then
              do ix=1,n_x_heg-1
                y=ymi+hh*ix
                p=k_fermi+alf_k_heg*sinh(y-x_fermi)
                xm(n_x_heg-ix)=(k2q2-p*p)/kq2
              enddo
              xm(0)=-1.d0
              xm(n_x_heg)=1.d0
              xw(0)=0.5d0*(xm(1)-xm(0))
              xw(n_x_heg)=0.5d0*(xm(n_x_heg)-xm(n_x_heg-1))
              do ix=1,n_x_heg-1
                xw(ix)=0.5d0*(xm(ix+1)-xm(ix-1))
              enddo
            endif
c -----------------------------------------------------------
            ax=0.d0
            ax0=0.d0
            do ix=0,n_x_heg
              x=xm(ix)
              p2=k2q2-kq2*x
              p=sqrt(p2)
              if(p>k_infinum) then
                do i_tau=0,n_tau
                  g(i_tau)=0.d0
                enddo
              else
                x=x_fermi+asinh((p-k_fermi)/alf_k_heg)
                i=x/step_k_x
                if(i*step_k_x<x) i=i+1
                if(i==0) i=1
                x=p-k_m_heg(i-1)
                do i_tau=0,n_tau
                  g(i_tau)=g2_spl(i,1,i_tau)+x*(g2_spl(i,2,i_tau)+x*(
     &                     g2_spl(i,3,i_tau)+x*g2_spl(i,4,i_tau)))
                enddo
              endif
              ax=ax+xw(ix)*g
              ax0=ax0+xw(ix)*g(0)
            enddo
            x=q2*q_w_heg(iq)
            do i_tau=0,n_tau
              i0_tau=min(i_tau,n_tau-i_tau)
              vrt_q0_nu0_heg(ik,i_tau,ispin)=
     &          vrt_q0_nu0_heg(ik,i_tau,ispin)+x*ax(i_tau)
     &                                        *wqt(iq,i0_tau)
            enddo
            vrt_stat_q0_nu0_heg(ik,ispin)=
     &        vrt_stat_q0_nu0_heg(ik,ispin)+8.d0*pi*ax0*q_w_heg(iq)
          enddo   !! over iq
	  enddo
      enddo
      vrt_q0_nu0_heg=vrt_q0_nu0_heg/4.d0/pi/pi
      vrt_stat_q0_nu0_heg=vrt_stat_q0_nu0_heg/4.d0/pi/pi
      do ispin=1,nspin
        do ik=0,n_k_heg
c          call from_tau_to_omega(vrt_om_q0_nu0_heg(ik,:,ispin),
c     &                           vrt_q0_nu0_heg(ik,:,ispin),1)
        enddo  !! over ik
      enddo
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_VRT_fine_OM_K_Re.'//ubi)
        do i_omega=0,n_omega
          aa=1.d0+real(vrt_om_q0_nu0_heg(0:4,i_omega,1))
     &           +vrt_stat_q0_nu0_heg(0:4,1)
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                                 (aa(ik),ik=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_fine_K_OM_Re.'//ubi)
        do ik=0,n_k_heg
          aa=1.d0+real(vrt_om_q0_nu0_heg(ik,0:4,1))
     &           +vrt_stat_q0_nu0_heg(ik,1)
          write(3,'(f9.3,5(1x,f12.7))')k_m_heg(ik),
     &                                 (aa(i_omega),i_omega=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_fine_OM_K_Im.'//ubi)
        do i_omega=0,n_omega
          aa=imag(vrt_om_q0_nu0_heg(0:4,i_omega,1))
          write(3,'(f9.3,5(1x,f12.7))')w_omega(i_omega)*evolt/2,
     &                                 (aa(ik),ik=0,4)
        enddo
        close(3)
        open(3,file=allfile(1:i_len)//'_VRT_fine_K_OM_Im.'//ubi)
        do ik=0,n_k_heg
          aa=imag(vrt_om_q0_nu0_heg(ik,0:4,1))
          write(3,'(f9.3,5(1x,f12.7))')k_m_heg(ik),
     &                                 (aa(i_omega),i_omega=0,4)
        enddo
        close(3)
      endif
      deallocate(g,ax,xm,xw,vrt_q0_nu0_heg,g2_om,g2_tau,g2_spl,
     &           vrt_om_q0_nu0_heg,vrt_stat_q0_nu0_heg)
	end