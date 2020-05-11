      subroutine heg1_pol
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer :: ispin,iq,ix,i,ii,i_tau,i1_tau
      real*8 :: k,q,x,p,k2,g1,green0,p2,q2,p2q2,pq2,de,piro4,hh,pmi,pma,
     &          ymi,yma,y,dg,dgx,dg1
      real*8, allocatable :: gx(:),g(:),xm(:),xw(:)
      allocate(gx(0:n_tau/2))
      allocate(g(0:n_tau/2))
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      piro4=4.d0*pi*ro_el_gas
      p_q_tau_heg=0.d0
      dp_q_tau_heg=0.d0
	do ispin=1,nspin
        do iq=0,n_q_heg
          q=q_m_heg(iq)
          q2=q*q
          do ii=0,n_p_heg
            p=p_m_heg(ii)
            p2=p*p
            p2q2=p2+q2
            pq2=2.d0*p*q
            gx=0.d0
            dgx=0.d0
c -------- Adaptive X-mesh ----------------------------------
            pmi=abs(p-q)
            pma=p+q
            ymi=x_fermi+asinh((pmi-k_fermi)/alf_k_heg)
            yma=x_fermi+asinh((pma-k_fermi)/alf_k_heg)
            hh=(yma-ymi)/n_x_heg
            xm=0.d0
            xw=2.d0/(n_x_heg+1)
            if(p>1.d-8.and.q>1.d-8) then
              do ix=1,n_x_heg-1
                y=ymi+hh*ix
                k=k_fermi+alf_k_heg*sinh(y-x_fermi)
                xm(n_x_heg-ix)=(p2q2-k*k)/pq2
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
            do ix=0,n_x_heg
              x=xm(ix)
              k2=p2q2-pq2*x
              k=sqrt(k2)
              if(k>k_infinum) then
                if(ubi_0=='dft') de=k2+vxc_heg(ispin)-chem_pot
                if(ubi_0/='dft') de=k2-piro4/k2-chem_pot
                do i_tau=0,n_tau/2
                  g(i_tau)=green0(de,tau_mesh(n_tau-i_tau))
                enddo
                dg=-de*g(0)
              else
                x=x_fermi+asinh((k-k_fermi)/alf_k_heg)
                i=x/step_k_x
                if(i*step_k_x<x) i=i+1
                if(i==0) i=1
                x=k-k_m_heg(i-1)
                do i_tau=0,n_tau/2
                  i1_tau=n_tau-i_tau
                  g(i_tau)=g_k_tau_spl(i,1,i1_tau,ispin)+x*(
     &                     g_k_tau_spl(i,2,i1_tau,ispin)+x*(
     &                     g_k_tau_spl(i,3,i1_tau,ispin)+x*
     &                     g_k_tau_spl(i,4,i1_tau,ispin)))
                enddo
                dg=dg_k_tau_spl(i,1,3,ispin)+x*(
     &             dg_k_tau_spl(i,2,3,ispin)+x*(
     &             dg_k_tau_spl(i,3,3,ispin)+x*
     &             dg_k_tau_spl(i,4,3,ispin)))
              endif
              gx=gx+xw(ix)*g
              dgx=dgx+xw(ix)*dg
            enddo
            if(p>k_infinum) then
              if(ubi_0=='dft') de=p2+vxc_heg(ispin)-chem_pot
              if(ubi_0/='dft') de=p2-piro4/p2-chem_pot
              do i_tau=0,n_tau/2
                g1=green0(de,tau_mesh(i_tau))
                p_q_tau_heg(iq,i_tau)=p_q_tau_heg(iq,i_tau)
     &                               -p_w_heg(ii)*p2*g1*gx(i_tau)
              enddo
              g1=green0(de,0.d0)
              dg1=-de*green0(de,0.d0)
              dp_q_tau_heg(iq)=dp_q_tau_heg(iq)
     &                        -p_w_heg(ii)*p2*(dg1*gx(0)-g1*dgx)
            else
              x=x_fermi+asinh((p-k_fermi)/alf_k_heg)
              i=x/step_k_x
              if(i*step_k_x<x) i=i+1
              if(i==0) i=1
              x=p-k_m_heg(i-1)
              do i_tau=0,n_tau/2
                g1=g_k_tau_spl(i,1,i_tau,ispin)+x*(
     &             g_k_tau_spl(i,2,i_tau,ispin)+x*(
     &             g_k_tau_spl(i,3,i_tau,ispin)+x*
     &             g_k_tau_spl(i,4,i_tau,ispin)))
                p_q_tau_heg(iq,i_tau)=p_q_tau_heg(iq,i_tau)
     &                               -p_w_heg(ii)*p2*g1*gx(i_tau)
              enddo
              g1=g_k_tau_spl(i,1,0,ispin)+x*(
     &           g_k_tau_spl(i,2,0,ispin)+x*(
     &           g_k_tau_spl(i,3,0,ispin)+x*
     &           g_k_tau_spl(i,4,0,ispin)))
              dg1=dg_k_tau_spl(i,1,1,ispin)+x*(
     &            dg_k_tau_spl(i,2,1,ispin)+x*(
     &            dg_k_tau_spl(i,3,1,ispin)+x*
     &            dg_k_tau_spl(i,4,1,ispin)))
              dp_q_tau_heg(iq)=dp_q_tau_heg(iq)
     &                        -p_w_heg(ii)*p2*(dg1*gx(0)-g1*dgx)
            endif
          enddo
	  enddo
      enddo
      p_q_tau_heg=p_q_tau_heg/4.d0/pi/pi
      if(nspin==1) p_q_tau_heg=p_q_tau_heg*2.d0
      dp_q_tau_heg=dp_q_tau_heg/4.d0/pi/pi
      if(nspin==1) dp_q_tau_heg=dp_q_tau_heg*2.d0
      deallocate(gx,g,xm,xw)
      call tau_to_nu_spl(p_q_nu_heg,p_q_tau_heg,dp_q_tau_heg,
     &                   n_q_heg+1,dp_tau_to_nu)
      call output_heg_p_tau(p_q_tau_heg,0)
	end