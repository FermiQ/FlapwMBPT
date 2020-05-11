      subroutine nu2_mesh_spl(v,wnu,tnu,num_nu)
      use heg_mod
      use manager_mod
      use parallel_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: v
      integer, intent(out) :: num_nu
      real*8, intent(out) :: wnu(n_nu2_max)
      complex*16, intent(out) :: tnu(n_nu2_max,2,ndim3_tau)
      logical :: left_exa
      integer :: i0,i3,i1,i2,k,i,i_nu,i_tau,key,ind_tau,it,i1_nu_m,j,
     &           ngeo1,nexa,ngeo2,nasy2,nas,i4
      real*8 :: hh,pib,w0,w1,w2,w3,c1,w4,aa,bb,xm,ct,st,hx,fac
      complex*16 :: cc
      real*8, allocatable :: mesh(:),tmp(:,:,:)
      pib=pi/betta_t
      hh=pib+pib
      nas=n_nu_geom+n_nu_asy
      allocate(mesh(0:nas))
      i1_nu_m=idnint(v/hh)
      key=0
      if(mod(i1_nu_m,2)==1) key=1
      i0=i1_nu_m/2
      if(key==1) i0=i0+1
      i1=i1_nu_m-n_nu_exa
      i2=i1_nu_m+n_nu_exa-1
      i3=i2+mats_num_nu(n_nu_exa+n_nu_geom)-n_nu_exa
      i4=i2+n_nu_maz-n_nu_exa
c --------------------------------------------------------------
      ngeo2=n_nu_geom
      nasy2=n_nu_asy
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
      w0=i0*hh
      w1=i1*hh
      w2=i2*hh
      w3=i3*hh
      w4=i4*hh
      k=0
      wnu(1)=w0
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
            tnu(1,it,ind_tau)=tnu(1,it,ind_tau)+c1*dcmplx(ct,-st)
          enddo
        enddo   
        wnu(2)=w0+hh
c ------ Left Geometric interval --------------------------------
        hx=(w1-wnu(2))/hh
        call geom_factor_finder(nas-1,hx,fac)
        hx=hh
        wnu(nas+1)=w1
        do i=1,nas-2
          k=nas-i+1
          wnu(k)=wnu(k+1)-hx
          j=idnint(wnu(k)/hh)
          wnu(k)=j*hh
          hx=hx*fac
        enddo
        do i_nu=0,nas
          mesh(i_nu)=wnu(i_nu+1)-v
        enddo
        allocate(tmp(nas+1,0:n_tau,2))
        call tau_from_freq_spl(n_tau,tau_mesh,nas,mesh,tmp,0)
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=ndim3_tau*me_t+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            if(it==2.and.i_tau==n_tau/2) cycle
            do i_nu=1,nas+1
              tnu(i_nu,it,ind_tau)=tnu(i_nu,it,ind_tau)
     &           +dcmplx(tmp(i_nu,i_tau,1),-tmp(i_nu,i_tau,2))
            enddo
          enddo
        enddo
        deallocate(tmp)
        k=nas
      endif
c ------ Exact interval --------------------------------
      do i=i1,i2
        k=k+1
        wnu(k)=hh*i
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=ndim3_tau*me_t+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            c1=(wnu(k)-v)*tau_mesh(i_tau)
            cc=dcmplx(cos(c1),sin(c1))
            if(abs(wnu(k)-v/2)<1.d-12) cc=cc/2
            tnu(k,it,ind_tau)=tnu(k,it,ind_tau)+cc
          enddo
        enddo
      enddo
c ------ Right geometric interval mesh --------------------------------
      hx=(w3-w2)/hh/2
      call geom_factor_finder(ngeo2,hx,fac)
      hx=hh*2
      do i=1,ngeo2
        k=k+1
        wnu(k)=wnu(k-1)+hx
        j=idnint(wnu(k)/hh)
        wnu(k)=j*hh
        hx=hx*fac
      enddo
c ------ Right asymptotic interval mesh -------------------------------
      wnu(k+nasy2+1)=w4
      wnu(k+nasy2)=w4-nu_last_step*hh
      bb=1.d0/w3**2
      aa=1.d0/w4**2-bb
      xm=(1.d0/wnu(k+nasy2)**2-bb)/aa
      hx=xm/nasy2
      do i=1,nasy2-1
        k=k+1
        xm=i*hx
        wnu(k)=1.d0/(aa*xm+bb)**(1.d0/2.d0)
        j=idnint(wnu(k)/hh)
        wnu(k)=j*hh
      enddo
      k=k+2
c ------ Right geometric interval weights -----------------------------
      do i_nu=0,ngeo2+2
        j=ngeo1+nexa-1+i_nu
        mesh(i_nu)=wnu(j)-v
      enddo
      allocate(tmp(ngeo2+3,0:n_tau,2))
      call tau_from_freq_spl(n_tau,tau_mesh,ngeo2+2,mesh,tmp,0)
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=ndim3_tau*me_t+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          if(it==2.and.i_tau==n_tau/2) cycle
          do i_nu=0,ngeo2+2
            j=ngeo1+nexa-1+i_nu
            tnu(j,it,ind_tau)=tnu(j,it,ind_tau)
     &         +dcmplx(tmp(i_nu+1,i_tau,1),-tmp(i_nu+1,i_tau,2))
          enddo
        enddo
      enddo
      deallocate(tmp)
c ------ Right asymptotic interval weights -----------------------------
      do i_nu=0,nasy2+1
        j=ngeo1+nexa+ngeo2+i_nu
        mesh(i_nu)=wnu(j)-v
      enddo
      allocate(tmp(nasy2+2,0:n_tau,2))
      call tau_from_freq_spl(n_tau,tau_mesh,nasy2+1,mesh,tmp,2)
      do ind_tau=1,ndim3_tau
        do it=1,2
          i_tau=ndim3_tau*me_t+ind_tau-1
          if(it==2) i_tau=n_tau-i_tau
          if(it==2.and.i_tau==n_tau/2) cycle
          do i_nu=0,nasy2+1
            j=ngeo1+nexa+ngeo2+i_nu
            tnu(j,it,ind_tau)=tnu(j,it,ind_tau)
     &         +dcmplx(tmp(i_nu+1,i_tau,1),-tmp(i_nu+1,i_tau,2))
          enddo
        enddo
      enddo
      deallocate(tmp,mesh)
c --------------------------------------------------------------
      num_nu=ngeo1+nexa+ngeo2+nasy2+1
      end
