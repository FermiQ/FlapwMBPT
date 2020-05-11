      subroutine self_energy_1_dyn_a(ispin,iq,ind_nu,nom,g_om_nu,
     &                               v2_q,s2,w2_q,vx_dyn,aw,z_red,a_red,
     &                               fif,phf,psi2)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      character*3, intent(in) :: psi2
      integer, intent(in) :: iq,ispin,ind_nu,nom,s2
	  real*8, intent(in) :: fif(maxel_red**2,n_pbmtm_red,nsort)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red),
     &                          w2_q(n_pbtot_red,n_pbtot_red),
     &                          vx_dyn(nrs_red,nrs_red,0:n_tau,2,
     &                                 nqdiv_c),
     &                          g_om_nu(nbndf_bnd,nbndf_bnd,2,nom,
     &                                  npnt_c),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
	  complex*16, intent(inout) :: aw(nra_red,nra_red,nqdiv_c,0:n_tau)
      integer :: k,kq,n,k0,nk,kq0,nkq,i,nn,i_omega,in,i_tau,mx2,iatom,
     &           isort,ind2,indp,ind,kf0,kfq0,iqa,na,nb0,nbas,ka,i0,ii,
     &           irl,ir,j
      real*8 :: v(3),dd1,dd2,dd3
      complex*16 :: c1
      complex*16, allocatable :: aa(:,:),tmp1(:,:,:),tmp2(:,:,:),
     &                           bb(:,:),vr(:,:,:,:),t0(:),
     &                           vrt(:,:,:),cc(:,:),w_mt(:),w_int(:),
     &                           t(:,:),t1(:,:),vw(:),gg(:,:)
      mx2=maxel_red**2
      iqa=k_a_from_c(iq)
      allocate(w_mt(n_lem2_red))
      allocate(w_int(nr_full_red))
      allocate(t0(nr_full_red))
      allocate(t(nrs_red,nra_red))
      allocate(t1(nra_red,nra_red))
      allocate(vw(n_pbtot_red))
      vw=v2_q(s2,:)+w2_q(s2,:)
	  do iatom=1,natom
	    isort=is(iatom)
	    n=lfunm_red(isort)
	    nn=n_pbmt0_red(isort)
	    ind2=io_lem2_red(iatom)
	    indp=iopb_red(iatom)
        call dgemm('n','t',2,n*n,nn,1.d0,vw(indp),2,fif(1,1,isort),mx2,
     &             0.d0,w_mt(ind2),2)
      enddo
      w_int=(0.d0,0.d0)
      do i=1,nplwgw_red(iqa)
        i0=indpw_gw_red(i,iq)
        call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
        w_int(ind)=vw(n_pbmt_red+i)
      enddo
      call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,w_int,-1)
      do ir=1,nr_full_red
        c1=conjg(phf(ir,iq))/sqrt(amega)
        w_int(ir)=c1*w_int(ir)
      enddo
      do kq=1,nqdiv_c  ! K+Q here
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
	    v=pnt_c(:,kq)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,k)
	    k=index_k1_c(k)     !! for K here
	    k0=i_kref_c(k)
        kf0=k_a_from_c(k0)
        ka=k_a_from_c(k)
	    nk=n_low_bnd(kf0,ispin)
        nb0=nbask_red(kf0)
        nbas=nb0/nrel
        na=nfun_red+nb0
c -------- VRT(tau) ---->  VRT(Omega) --------------------------
        allocate(vr(nkq,na,0:n_tau,2))
        do in=1,2 
          do i_tau=0,n_tau
c --- Transform VRT_(tL)r_(t'L')r' to VRT_lambda_(tL')G' -----------
            call vrt_r_r1_to_bnd_g1(kq,k,nkq,na,z_red(1,1,kq),
     &                              a_red(1,1,kq),phf,
     &                              vx_dyn(1,1,i_tau,in,kq),
     &                              vr(1,1,i_tau,in))
          enddo
        enddo
        allocate(vrt(nkq,na,nom))
        allocate(aa(nkq,na))	 
        allocate(bb(nkq,nrs_red))
        allocate(cc(nrs_red,nkq))
        allocate(tmp1(nrs_red,na,n_omega1_max))
        allocate(tmp2(nrs_red,na,0:n_tau))
        allocate(gg(nkq,nkq))
        do in=1,2 
          call from_tau_nu_to_omega_nu_d(na*nkq,ind_nu,nom,vrt,vr,in)
c ------------------------------------------------------------------          
	      tmp1=(0.d0,0.d0)
	      do i_omega=1,nom
	        aa=vrt(:,:,i_omega)
	        if(in==1) then  ! f(-w+v) here
	          do j=1,nkq
	            do i=1,nkq
	              gg(i,j)=conjg(g_om_nu(j,i,2,i_omega,kq0))
	            enddo
	          enddo
	        else
	          gg=g_om_nu(1:nkq,1:nkq,1,i_omega,kq0)
	        endif
c --- Transform G_bnd_bnd' to G_(tL)r_bnd' -------------------------
            call g_bnd_bnd1_to_r_bnd1(kq,nkq,z_red(1,1,kq),
     &                                a_red(1,1,kq),phf,gg,cc)
            call zgemm('n','n',nrs_red,na,nkq,(1.d0,0.d0),cc,nrs_red,
     &                 aa,nkq,(0.d0,0.d0),tmp1(1,1,i_omega),nrs_red)
          enddo
          call from_omega_nu_to_tau_nu_aa(nrs_red*na,ind_nu,tmp1,tmp2,
     &                                    in)
c -------------------------------------------------------------------
          dd3=wgt_c(iq)*tau_from_nu_d(ind_nu,0,1)
          do i_tau=0,n_tau
	        do iatom=1,natom
	          isort=is(iatom)
	          n=lfunm_red(isort)
	          ind2=io_lem2_red(iatom)
	          ind=io_lem_red(iatom)
	          call zgemm('t','n',n,na,n,(1.d0,0.d0),w_mt(ind2),n,
     &	                 tmp2(ind,1,i_tau),nrs_red,(0.d0,0.d0),t(ind,1),
     &                   nrs_red)
            enddo
            do i=1,nr_full_red
              t(nfun_red+i,1:na)=tmp2(nfun_red+i,1:na,i_tau)*w_int(i)
            enddo
            t1(1:nfun_red,1:na)=t(1:nfun_red,1:na)
c --- Transform SIG_r_G' to SIG_G_G' -------------------------------
            do j=1,na
              do irl=1,nrel
                ind=nfun_red+(irl-1)*nr_full_red
                do ir=1,nr_full_red
                  c1=conjg(phf(ir,k))
                  t0(ir)=c1*t(ind+ir,j)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),1,t0,
     &                    -1)
                ii=nfun_red+(irl-1)*nbas
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  t1(ii+i,j)=t0(ind)/nr_full_red*sqrt(amega)
                enddo
              enddo   !! over irl
            enddo
c ------------------------------------------------------------------
            dd1=wgt_c(iq)*tau_from_nu_d(ind_nu,i_tau,1)
            dd2=wgt_c(iq)*tau_from_nu_d(ind_nu,i_tau,2)
            if(in==1) then
              aw(1:na,1:na,k,i_tau)=aw(1:na,1:na,k,i_tau)
     &                               +dcmplx(dd1,dd2)*t1(1:na,1:na)
            else
              aw(1:na,1:na,k,i_tau)=aw(1:na,1:na,k,i_tau)
     &                               +dd3*t1(1:na,1:na)
            endif
          enddo  !! over i_tau
        enddo   !! over in
        deallocate(tmp1,tmp2,vr,vrt,aa,bb,cc,gg)
      enddo  !! over kq
      deallocate(w_mt,w_int,t,t1,vw,t0)
      end
