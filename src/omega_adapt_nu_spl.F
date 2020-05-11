      subroutine omega_adapt_nu_spl(v,wom,tom,num_om)
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: v
      real*8, intent(out) :: wom(n_omega1_max)
      complex*16, intent(out) :: tom(n_omega1_max,2,ndim3_tau)
      integer, intent(out) :: num_om
      logical :: left_exa
      integer :: i_omega,i0,i3,i1,i2,k,i,i_tau,key,i_nu_m,ind_tau,it,
     &           nas,j,i4,ngeo1,nexa,ngeo2,nasy2
      real*8 :: hh,pib,w0,w1,w2,w3,c1,ct,st,fac,hx,w4,aa,bb,xm
      complex*16 :: cc
      real*8, allocatable :: tmp(:,:,:),mesh(:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_omega_geom+n_omega_asy
      allocate(mesh(0:nas))
      i_nu_m=idnint(v/hh)
      i0=i_nu_m/2
      i1=i_nu_m-n_omega_exa
      i2=i_nu_m+n_omega_exa
      i3=i2+mats_num_omega(n_omega_exa+n_omega_geom)-n_omega_exa
      i4=i2+n_omega_maz-n_omega_exa
c --------------------------------------------------------------
      ngeo2=n_omega_geom
      nasy2=n_omega_asy
      if(i1-i0<=nas) then
        i1=i0
        ngeo1=0
        left_exa=.true.
      else
        ngeo1=nas
        left_exa=.false.
      endif
      nexa=i2-i1+1
c --------------------------------------------------------------
      w0=pib+i0*hh
      w1=pib+i1*hh
      w2=pib+i2*hh
      w3=pib+i3*hh
      w4=pib+i4*hh
      wom(1)=w0
      key=0
      if(mod(i_nu_m,2)==1) key=1
      k=0
      if(.not.left_exa) then
c ------------ The first point -------------------------------------
        c1=1.d0
        if(key==0) c1=0.5d0
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=ndim3_tau*me_t+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            if(it==2.and.i_tau==n_tau/2) cycle
            ct=cos((w0-v)*tau_mesh(i_tau))
            st=sin((w0-v)*tau_mesh(i_tau))
            tom(1,it,ind_tau)=tom(1,it,ind_tau)+c1*dcmplx(ct,-st)
          enddo
        enddo   
        wom(2)=w0+hh
c ------ Left Geometric interval --------------------------------
        hx=(w1-wom(2))/hh
        call geom_factor_finder(nas-1,hx,fac)
        hx=hh
        wom(nas+1)=w1
        do i=1,nas-2
          k=nas-i+1
          wom(k)=wom(k+1)-hx
          j=idnint((wom(k)-pib)/hh)
          wom(k)=(2*j+1)*pib
          hx=hx*fac
        enddo
        do i_omega=0,nas
          mesh(i_omega)=wom(i_omega+1)-v
        enddo
        allocate(tmp(nas+1,0:n_tau,2))
        call tau_from_freq_spl(n_tau,tau_mesh,nas,mesh,tmp,0)
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=ndim3_tau*me_t+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            if(it==2.and.i_tau==n_tau/2) cycle
            do i_omega=1,nas+1
              tom(i_omega,it,ind_tau)=tom(i_omega,it,ind_tau)
     &           +dcmplx(tmp(i_omega,i_tau,1),-tmp(i_omega,i_tau,2))
            enddo
          enddo
        enddo
        deallocate(tmp)
        k=nas
      endif
c ------ Exact interval --------------------------------
      do i=i1,i2
        k=k+1
        wom(k)=pib+hh*i
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=ndim3_tau*me_t+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            c1=-(wom(k)-v)*tau_mesh(i_tau)
            cc=dcmplx(cos(c1),sin(c1))
            if(abs(wom(k)-v/2)<1.d-12) cc=cc/2
            tom(k,it,ind_tau)=tom(k,it,ind_tau)+cc
          enddo
        enddo
      enddo
c ------ Right geometric interval mesh --------------------------------
      hx=(w3-w2)/hh/2
      call geom_factor_finder(ngeo2,hx,fac)
      hx=hh*2
      do i=1,ngeo2
        k=k+1
        wom(k)=wom(k-1)+hx
        j=idnint((wom(k)-pib)/hh)
        wom(k)=(2*j+1)*pib
        hx=hx*fac
      enddo
c ------ Right asymptotic interval mesh -------------------------------
      wom(k+nasy2+1)=w4
      wom(k+nasy2)=w4-omega_last_step*hh
      bb=1.d0/w3**2
      aa=1.d0/w4**2-bb
      xm=(1.d0/wom(k+nasy2)**2-bb)/aa
      hx=xm/nasy2
      do i=1,nasy2-1
        k=k+1
        xm=i*hx
        wom(k)=1.d0/(aa*xm+bb)**(1.d0/2.d0)
        j=idnint((wom(k)-pib)/hh)
        wom(k)=(2*j+1)*pib
      enddo
      k=k+2
c ------ Right geometric interval weights -----------------------------
      do i_omega=0,ngeo2+2
        j=ngeo1+nexa-1+i_omega
        mesh(i_omega)=wom(j)-v
      enddo
      allocate(tmp(ngeo2+3,0:n_tau,2))
      call tau_from_freq_spl(n_tau,tau_mesh,ngeo2+2,mesh,tmp,0)
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=ndim3_tau*me_t+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          if(it==2.and.i_tau==n_tau/2) cycle
          do i_omega=0,ngeo2+2
            j=ngeo1+nexa-1+i_omega
            tom(j,it,ind_tau)=tom(j,it,ind_tau)
     &         +dcmplx(tmp(i_omega+1,i_tau,1),-tmp(i_omega+1,i_tau,2))
          enddo
        enddo
      enddo
      deallocate(tmp)
c ------ Right asymptotic interval weights -----------------------------
      do i_omega=0,nasy2+1
        j=ngeo1+nexa+ngeo2+i_omega
        mesh(i_omega)=wom(j)-v
      enddo
      allocate(tmp(nasy2+2,0:n_tau,2))
      call tau_from_freq_spl(n_tau,tau_mesh,nasy2+1,mesh,tmp,2)
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=ndim3_tau*me_t+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          if(it==2.and.i_tau==n_tau/2) cycle
          do i_omega=0,nasy2+1
            j=ngeo1+nexa+ngeo2+i_omega
            tom(j,it,ind_tau)=tom(j,it,ind_tau)
     &         +dcmplx(tmp(i_omega+1,i_tau,1),-tmp(i_omega+1,i_tau,2))
          enddo
        enddo
      enddo
      deallocate(tmp,mesh)
c --------------------------------------------------------------
      num_om=ngeo1+nexa+ngeo2+nasy2+1
      end
