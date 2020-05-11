	real*8 function heg_sig_x_exact(k)
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	real*8, intent(in) :: k
	real*8 :: y,s,pi2a
      pi2a=(pi+pi)/par
	y=k/k_fermi
	if(k<1.d-9) then
	  s=-2.d0
	else if(abs(y-1.d0)<1.d-12) then
	  s=-1.d0
	else
	  s=abs((1.d0+y)/(1.d0-y))
	  s=(1.d0-y*y)*log(s)/2.d0/y
	  s=-(1.d0+s)
	endif
	heg_sig_x_exact=2.d0*k_fermi*s/pi
      end
      
      
      subroutine heg1_sig_x
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer :: ik,ispin,iq,ix,i
      real*8 :: ax,k,q,k2q2,kq2,x,p,k2,g,fermi_dirac,p2,piro4,hh,pmi,
     &          pma,ymi,yma,y
      real*8, allocatable :: xm(:),xw(:)
      piro4=4.d0*pi*ro_el_gas
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      sigx_k_heg=0.d0
	do ispin=1,nspin
	  do ik=0,n_k_heg
          k=k_m_heg(ik)
          k2=k*k
          do iq=0,n_p_heg
            q=p_m_heg(iq)
            k2q2=k2+q*q
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
            do ix=0,n_x_heg
              x=xm(ix)
              p2=k2q2-kq2*x
              if(p2<0.d0) p2=0.d0
              p=sqrt(p2)
              if(p>k_infinum) then
                if(ubi_0=='dft') g=-fermi_dirac(p2+vxc_heg(ispin)
     &                                         -chem_pot)
                if(ubi_0/='dft') g=-fermi_dirac(p2-piro4/p2-chem_pot)
              else
                x=x_fermi+asinh((p-k_fermi)/alf_k_heg)
                i=x/step_k_x
                if(i*step_k_x<x) i=i+1
                if(i==0) i=1
                i=min(i,n_k_heg)
                x=p-k_m_heg(i-1)
                g=g_k_tau_spl(i,1,n_tau,ispin)+x*(
     &            g_k_tau_spl(i,2,n_tau,ispin)+x*(
     &            g_k_tau_spl(i,3,n_tau,ispin)+x*
     &            g_k_tau_spl(i,4,n_tau,ispin)))
              endif
              ax=ax+xw(ix)*g
            enddo
            sigx_k_heg(ik,ispin)=sigx_k_heg(ik,ispin)+p_w_heg(iq)*ax
          enddo
	  enddo
      enddo
      sigx_k_heg=sigx_k_heg*2.d0/pi
      deallocate(xm,xw)
	end