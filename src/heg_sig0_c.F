	subroutine heg_sig0_c
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
	integer :: ik,ispin,iq,ix,i_omega
      real*8 :: k,q,k2q2,kq2,x,p,k2,p2,q2,piro4,hh,pmi,pma,ymi,yma,y,
     &          epm,kap,aa,bb,ekap1,ekap2,w
      real*8, allocatable :: xm(:),xw(:)
      complex*16, allocatable :: ax(:),g(:)
      piro4=4.d0*pi*ro_el_gas
      allocate(g(0:n_omega))
      allocate(ax(0:n_omega))
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      sigc0_k_omega_heg=(0.d0,0.d0)
	do ispin=1,nspin
	  do ik=0,n_k_heg
          k=k_m_heg(ik)
          k2=k*k
          do iq=0,n_q_heg
            q=q_m_heg(iq)
            q2=q*q
            k2q2=k2+q2
            kq2=2.d0*k*q
            kap=sqrt(p_as_nu_par(iq,1)**2+16.d0*pi*ro_el_gas)
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
            ax=(0.d0,0.d0)
            do ix=0,n_x_heg
              x=xm(ix)
              p2=k2q2-kq2*x
              p=sqrt(p2)
              epm=p2+vxc_heg(ispin)-chem_pot_lda
              call sig0_ab(epm,kap,aa,bb)
              ekap1=-epm-kap
              ekap2=kap-epm
              do i_omega=0,n_omega
                w=w_omega(i_omega)
                g(i_omega)=aa/dcmplx(ekap1,w)+bb/dcmplx(ekap2,w)
              enddo
              ax=ax+xw(ix)*g
            enddo
            x=q_w_heg(iq)/kap
            do i_omega=0,n_omega
              sigc0_k_omega_heg(ik,i_omega,ispin)=
     &          sigc0_k_omega_heg(ik,i_omega,ispin)+x*ax(i_omega)
            enddo
          enddo   !! over iq
	  enddo
      enddo
      sigc0_k_omega_heg=sigc0_k_omega_heg*16.d0*ro_el_gas
      deallocate(g,ax,xm,xw)
      end
      
      
      subroutine sig0_ab(epm,kap,aa,bb)
	use manager_mod
      implicit none
	real*8, intent(in) :: kap,epm
      real*8, intent(out) :: aa,bb
	real*8 :: a,b,t1,t2,t3,c
      if(epm<0.d0) then
	  t1=-epm*betta_t
	  t2=kap*betta_t
	  t3=(kap-epm)*betta_t
	  a=0.d0
	  if(t1<100.d0) a=exp(-t1)
	  if(t2<100.d0) a=a+exp(-t2)
	  b=1.d0
	  if(t3<100.d0) b=b+exp(-t3)
        c=1.d0
	  if(t1<100.d0) c=c+exp(-t1)
	  if(t2<100.d0) c=c-exp(-t2)
	  if(t3<100.d0) c=c-exp(-t3)
      else
	  t1=epm*betta_t
	  t2=kap*betta_t
	  t3=(kap+epm)*betta_t
	  a=1.d0
	  if(t3<100.d0) a=a+exp(-t3)
	  b=0.d0
	  if(t1<100.d0) b=b+exp(-t1)
	  if(t2<100.d0) b=b+exp(-t2)
        c=1.d0
	  if(t1<100.d0) c=c+exp(-t1)
	  if(t2<100.d0) c=c-exp(-t2)
	  if(t3<100.d0) c=c-exp(-t3)
      endif
      aa=a/c
      bb=b/c
      end
      
      
      subroutine heg_sig0_tau(wqt)
	use heg_mod
	use manager_mod
	use solid_mod
	use units_mod
	implicit none
      real*8, intent(in) :: wqt(0:n_q_heg,0:n_tau/2)
	integer :: ik,ispin,iq,ix,i_tau,i0_tau
      real*8 :: k,q,k2q2,kq2,x,p,k2,green0,p2,de,q2,piro4,hh,pmi,pma,
     &          ymi,yma,y
      real*8, allocatable :: g(:),ax(:),xm(:),xw(:)
      piro4=4.d0*pi*ro_el_gas
      allocate(g(0:n_tau))
      allocate(ax(0:n_tau))
      allocate(xm(0:n_x_heg))
      allocate(xw(0:n_x_heg))
      sigc0_k_tau_heg=0.d0
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
            do ix=0,n_x_heg
              x=xm(ix)
              p2=k2q2-kq2*x
              p=sqrt(p2)
              de=p2+vxc_heg(ispin)-chem_pot_lda
              do i_tau=0,n_tau
                g(i_tau)=green0(de,tau_mesh(i_tau))
              enddo
              ax=ax+xw(ix)*g
            enddo
            x=q2*q_w_heg(iq)
            do i_tau=0,n_tau
              i0_tau=min(i_tau,n_tau-i_tau)
              sigc0_k_tau_heg(ik,i_tau,ispin)=
     &          sigc0_k_tau_heg(ik,i_tau,ispin)-x*ax(i_tau)
     &                                        *wqt(iq,i0_tau)
            enddo
          enddo   !! over iq
	  enddo
      enddo
      sigc0_k_tau_heg=sigc0_k_tau_heg/4.d0/pi/pi
      deallocate(g,ax,xm,xw)
	end