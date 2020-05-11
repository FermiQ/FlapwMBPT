	subroutine heg1_sig_c(wqt)
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
      real*8, intent(in) :: wqt(0:n_q_heg,0:n_tau/2)
	integer :: ik,ispin,iq,ix,i,i_tau,i0_tau
      real*8 :: k,q,k2q2,kq2,x,p,k2,green0,p2,de,q2,piro4,hh,pmi,pma,
     &          ymi,yma,y,dg(3),dax(3)
      real*8, allocatable :: g(:),ax(:),xm(:),xw(:)
      piro4=4.d0*pi*ro_el_gas
      allocate(g(0:n_tau))
      allocate(ax(0:n_tau))
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      sigc_k_tau_heg=0.d0
      dsigc_k_tau_heg=0.d0
	do ispin=1,nspin
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
            dax=0.d0
            do ix=0,n_x_heg
              x=xm(ix)
              p2=k2q2-kq2*x
              p=sqrt(p2)
              if(p>k_infinum) then
                if(ubi_0=='dft') de=p2+vxc_heg(ispin)-chem_pot
                if(ubi_0/='dft') de=p2-chem_pot-piro4/p2
                do i_tau=0,n_tau
                  g(i_tau)=green0(de,tau_mesh(i_tau))
                enddo
                dg(1)=-de*g(0)
                dg(2)=-de*g(n_tau/2)
                dg(3)=-de*g(n_tau)
              else
                x=x_fermi+asinh((p-k_fermi)/alf_k_heg)
                i=x/step_k_x
                if(i*step_k_x<x) i=i+1
                if(i==0) i=1
                i=min(i,n_k_heg)
                x=p-k_m_heg(i-1)
                do i_tau=0,n_tau
                  g(i_tau)=g_k_tau_spl(i,1,i_tau,ispin)+x*(
     &                     g_k_tau_spl(i,2,i_tau,ispin)+x*(
     &                     g_k_tau_spl(i,3,i_tau,ispin)+x*
     &                     g_k_tau_spl(i,4,i_tau,ispin)))
                enddo
                do i_tau=1,3
                  dg(i_tau)=dg_k_tau_spl(i,1,i_tau,ispin)+x*(
     &                      dg_k_tau_spl(i,2,i_tau,ispin)+x*(
     &                      dg_k_tau_spl(i,3,i_tau,ispin)+x*
     &                      dg_k_tau_spl(i,4,i_tau,ispin)))
                enddo
              endif
              ax=ax+xw(ix)*g
              dax=dax+xw(ix)*dg
            enddo
            x=q2*q_w_heg(iq)
            do i_tau=0,n_tau
              i0_tau=min(i_tau,n_tau-i_tau)
              sigc_k_tau_heg(ik,i_tau,ispin)=
     &          sigc_k_tau_heg(ik,i_tau,ispin)-x*ax(i_tau)
     &                                        *wqt(iq,i0_tau)
            enddo
            dsigc_k_tau_heg(ik,1,ispin)=
     &          dsigc_k_tau_heg(ik,1,ispin)
     &         -x*(dax(1)*wqt(iq,0)+ax(0)*dw_q_tau_heg(iq))
            dsigc_k_tau_heg(ik,2,ispin)=
     &          dsigc_k_tau_heg(ik,2,ispin)
     &         -x*dax(2)*wqt(iq,n_tau/2)
            dsigc_k_tau_heg(ik,3,ispin)=
     &          dsigc_k_tau_heg(ik,3,ispin)
     &         -x*(dax(3)*wqt(iq,0)+ax(n_tau)*dw_q_tau_heg(iq))
          enddo   !! over iq
	  enddo
      enddo
      sigc_k_tau_heg=sigc_k_tau_heg/4.d0/pi/pi
      dsigc_k_tau_heg=dsigc_k_tau_heg/4.d0/pi/pi
      deallocate(g,ax,xm,xw)
      end
      
      
      subroutine heg_sigc_gw_cube
	use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      use vertex_mod
	implicit none
	integer :: i_tau,i0_tau,ispin,ik,iq,iq0,kq,kq0
      real*8 :: v(3)
	sigc_k_tau_heg_c=0.d0
	do i_tau=0,n_tau
	  i0_tau=min(i_tau,n_tau-i_tau)
	  do ispin=1,nspin
          do ik=1,nstar_c     
	      do iq=1,nrr_div_c
	        v=q_irr_c(:,ik)-q_sph_c(:,iq)
	        call zone1_number(v,rb0_c,nr_div_c,kq)
	        iq0=i_kref_c(iq)
	        kq0=i_kref_c(kq)
	        sigc_k_tau_heg_c(i_tau,ik,ispin)=
     &          sigc_k_tau_heg_c(i_tau,ik,ispin)
     &         -g_k_tau_heg_c(i_tau,kq0,ispin)*w_q_tau_heg_c(i0_tau,iq0)
	      enddo
	    enddo
	  enddo
      enddo
      sigc_k_tau_heg_c=sigc_k_tau_heg_c/amega/nqdiv_c
	end