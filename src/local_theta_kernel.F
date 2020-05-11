      subroutine local_theta_kernel     
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: mx2,iatom,ispin,isort,ii,i,j,ij,n,np,ind_tau,i_tau,
     &           ind_omega,i_omega,it,om1,om2,nu1,nu2
      real*8, allocatable :: g_tau(:,:,:,:,:,:),fif(:,:,:,:),
     &                       vloc(:,:,:),wloc(:,:,:,:),gt(:,:,:,:),
     &                       v2(:,:),w2(:,:,:),v0(:,:,:,:)
      complex*16, allocatable :: g_omega(:,:,:,:,:),z_red(:,:,:,:),
     &                           go(:,:,:,:),th(:,:,:,:,:,:,:)
      mx2=maxel_red**2
      allocate(fif(mx2,n_pbmtm_red,nsort,nspin))
      allocate(z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin))
      call aux_vrt_loc(fif,z_red)
      allocate(g_tau(maxel_red,maxel_red,2,ndim3_tau,nspin,natom))
      allocate(g_omega(maxel_red,maxel_red,ndim3_omega,nspin,natom))
      call g_loc_from_full(g_omega,g_tau,z_red)
      allocate(vloc(n_pbmtm_red,n_pbmtm_red,natom))
      allocate(wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu,natom))
      call w_loc_from_red(vloc,wloc)
      om1=-1
      om2=0
      nu1=0
      nu2=10
      do iatom=1,natom
        isort=is(iatom)
        n=lfunm_red(isort)
        if(n==0) cycle
        np=n_pbmt0_red(isort)
        allocate(gt(n,n,0:n_tau,nspin))
        gt=0.d0
        do ind_tau=1,ndim3_tau
          do it=1,2
            i_tau=me_t*ndim3_tau+ind_tau-1
            if(it==2) i_tau=n_tau-i_tau
            do ispin=1,nspin
              gt(:,:,i_tau,ispin)=g_tau(1:n,1:n,it,ind_tau,ispin,iatom)
            enddo
          enddo
        enddo
        allocate(go(n,n,0:n_omega,nspin))
        go=(0.d0,0.d0)
        do ind_omega=1,ndim3_omega
          i_omega=me_t*ndim3_omega+ind_omega-1
          do ispin=1,nspin
            go(:,:,i_omega,ispin)=g_omega(1:n,1:n,ind_omega,ispin,iatom)
          enddo
        enddo
	  if(nproc_t/=1) then
          call dgop(gt,n*n*(n_tau+1)*nspin,'  +',comm_t)
          call dgop(go,2*n*n*(n_omega+1)*nspin,'  +',comm_t)
        endif
        allocate(w2(np,np,0:n_nu))
        allocate(v2(np,np))
        v2=vloc(1:np,1:np,iatom)
        w2=wloc(1:np,1:np,:,iatom)
        allocate(v0(n,n,np,nspin))
        do ispin=1,nspin
          do ii=1,np
            ij=0
            do j=1,n
              do i=1,n
                ij=ij+1
                v0(i,j,ii,ispin)=fif(ij,ii,isort,ispin)
              enddo
            enddo
          enddo
        enddo
        allocate(th(2*n,2*n,2*n,2*n,om1:om2,om1:om2,nu1:nu2))
        call theta_kernel_atom(n,np,gt,go,v2,w2,v0,om1,om2,nu1,nu2,th)
        deallocate(gt,go,v2,w2,v0,th)
      enddo
      deallocate(g_tau,g_omega,fif,z_red,vloc,wloc)
      call ending
      end
      
      
      subroutine theta_kernel_atom(n,np,gt,go,v2,w2,v0,om1,om2,nu1,nu2,
     &                             th)     
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,np,om1,om2,nu1,nu2
      real*8, intent(in) :: gt(n,n,0:n_tau,nspin),w2(np,np,0:n_nu),
     &                      v2(np,np),v0(n,n,np,nspin)
      complex*16, intent(in) :: go(n,n,0:n_omega,nspin)
      complex*16, intent(out) :: th(2*n,2*n,2*n,2*n,om1:om2,om1:om2,
     &                              nu1:nu2) 
      integer :: ispin,jspin,i_nu,iom1,iom2,nn,n2,ij,j,i,it
      real*8 :: nu,o1,o2,pib
      real*8, allocatable :: v4(:,:,:,:),w4(:,:,:,:,:),tmp(:,:),
     &                       wspl_geom(:,:,:,:,:),wspl_asy(:,:,:,:,:),
     &                       trn(:),g0(:,:),gspl_geom(:,:,:,:,:,:),
     &                       gasy(:,:,:,:),gspl_asy(:,:,:,:,:,:),
     &                       gexa(:,:,:,:,:),wexa(:,:,:,:),wasy(:,:,:)
      complex*16, allocatable :: th0(:,:,:,:),th1(:,:,:,:,:,:,:)
      nn=n*n
      n2=nn*nn
      allocate(v4(nn,nn,nspin,nspin))
      allocate(w4(nn,nn,0:n_nu,nspin,nspin))
      allocate(tmp(nn,np))
      do jspin=1,nspin
        do ispin=1,nspin
	    call dgemm('n','n',nn,np,np,1.d0,v0(1,1,1,ispin),nn,v2,np,
     &               0.d0,tmp,nn)
	    call dgemm('n','t',nn,nn,np,1.d0,tmp,nn,v0(1,1,1,jspin),
     &               nn,0.d0,v4(1,1,ispin,jspin),nn)
          do i_nu=0,n_nu
            call dgemm('n','n',nn,np,np,1.d0,v0(1,1,1,ispin),nn,
     &                 w2(1,1,i_nu),np,0.d0,tmp,nn)
	      call dgemm('n','t',nn,nn,np,1.d0,tmp,nn,v0(1,1,1,jspin),
     &                 nn,0.d0,w4(1,1,i_nu,ispin,jspin),nn)
          enddo
        enddo
      enddo
      deallocate(tmp)
c ---------- Prepare interpolation --------------------------------
      allocate(gexa(0:n_omega_exa,2,n,n,nspin))
      allocate(gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin))
      allocate(gspl_asy(n_omega_asy+1,4,2,n,n,nspin))
      allocate(gasy(2,n,n,nspin))
      allocate(g0(0:n_omega,2))
      do ispin=1,nspin
        do j=1,n
          do i=1,n
            g0(:,1)=real(go(i,j,:,ispin))
            g0(:,2)=imag(go(i,j,:,ispin))
            do it=1,2
              call spline_inhmg(w_omega(n_omega_exa-1),
     &                          g0(n_omega_exa-1,it),
     &                          gspl_geom(0,1,it,i,j,ispin),
     &                          gspl_geom(0,2,it,i,j,ispin),
     &                          gspl_geom(0,3,it,i,j,ispin),
     &                          gspl_geom(0,4,it,i,j,ispin),
     &                          n_omega_geom+3,0,0.d0,0.d0)
              call spline_inhmg(xm_omega(0),
     &                          g0(n_omega_exa+n_omega_geom,it),
     &                          gspl_asy(1,1,it,i,j,ispin),
     &                          gspl_asy(1,2,it,i,j,ispin),
     &                          gspl_asy(1,3,it,i,j,ispin),
     &                          gspl_asy(1,4,it,i,j,ispin),
     &                          n_omega_asy+1,0,0.d0,0.d0)
              gasy(it,i,j,ispin)=g0(n_omega,it)
              gexa(:,it,i,j,ispin)=g0(0:n_omega_exa,it)
            enddo
          enddo
        enddo
      enddo
      deallocate(g0)
c -----------------------------------------------------------------
      allocate(wexa(0:n_nu_exa,n2,nspin,nspin))
      allocate(wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin))
      allocate(wspl_asy(n_nu_asy+1,4,n2,nspin,nspin))
      allocate(wasy(n2,nspin,nspin))
      do jspin=1,nspin
        do ispin=1,nspin
          ij=0
          do j=1,nn
            do i=1,nn
              ij=ij+1
              call spline_inhmg(w_nu(n_nu_exa-1),
     &                          w4(i,j,n_nu_exa-1:n_nu_exa+n_nu_geom+2,
     &                             ispin,jspin),
     &                          wspl_geom(0,1,ij,ispin,jspin),
     &                          wspl_geom(0,2,ij,ispin,jspin),
     &                          wspl_geom(0,3,ij,ispin,jspin),
     &                          wspl_geom(0,4,ij,ispin,jspin),
     &                          n_nu_geom+3,0,0.d0,0.d0)
              call spline_inhmg(xm_nu(0),
     &                          w4(i,j,n_nu_exa+n_nu_geom:n_nu_exa
     &                              +n_nu_geom+n_nu_asy+1,ispin,jspin),
     &                          wspl_asy(1,1,ij,ispin,jspin),
     &                          wspl_asy(1,2,ij,ispin,jspin),
     &                          wspl_asy(1,3,ij,ispin,jspin),
     &                          wspl_asy(1,4,ij,ispin,jspin),
     &                          n_nu_asy+1,0,0.d0,0.d0)
              wexa(:,ij,ispin,jspin)=w4(i,j,0:n_nu_exa,ispin,jspin)
              wasy(ij,ispin,jspin)=w4(i,j,n_nu,ispin,jspin)
            enddo
          enddo
        enddo
      enddo
c ---------------------------------------------------------------
      pib=pi/betta_t
      allocate(trn(0:n_nu))
      trn=0.d0
      trn(0)=tau_from_nu(0,0,1)/2.d0
      do i_nu=1,n_nu
        trn(i_nu)=tau_from_nu(i_nu,0,1)/4.d0
      enddo
      th=(0.d0,0.d0)
      allocate(th0(2*n,2*n,2*n,2*n))
      allocate(th1(2*n,2*n,2*n,2*n,om1:om2,om1:om2,nu1:nu2))
      th1=(0.d0,0.d0)
      do i_nu=nu1,nu2
        nu=2*i_nu*pib
        do iom1=om1,om2
          o1=(2*iom1+1)*pib
          do iom2=om1,om2
            o2=(2*iom2+1)*pib
            call kernel_loc_1_1(n,n2,th0,nu,v4,wexa,wasy,wspl_geom,
     &                          wspl_asy)
            th1(:,:,:,:,iom1,iom2,i_nu)=th1(:,:,:,:,iom1,iom2,i_nu)+th0
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
            call kernel_loc_1_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
            call kernel_loc_2_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
            call kernel_loc_3_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
            call kernel_loc_4_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
            call kernel_loc_5_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
            th(:,:,:,:,iom1,iom2,i_nu)=th(:,:,:,:,iom1,iom2,i_nu)+th0
          enddo
        enddo
      enddo
      call output_theta(n,th1(1,1,1,1,0,0,0),'1_1')
      call output_theta(n,th(1,1,1,1,0,0,0),'1_2')
      deallocate(v4,w4,gspl_geom,gspl_asy,gasy,gexa,wspl_geom,
     &           wspl_asy,wexa,wasy,th0,trn,th1)
      end