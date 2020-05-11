      subroutine t_rr_space(qt_pw,vloc,ind_nu,nomnu,wloc)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nomnu
      real*8, intent(in) :: wloc(n_pbmtm_red,n_pbmtm_red,0:n_nu,natom),
     &                      vloc(n_pbmtm_red,n_pbmtm_red,natom)
      complex*16, intent(inout) :: qt_pw(n_pbmt_red,n_pbmt_red,
     &                                   nqdiv_c,2,2,0:n_tau/2)
      integer :: j,i,nn,i_nu1,in,i_tau,i_nu,i0_tau,it,il,iatom,isort,
     &           jatom,jsort,ir,n1,jnd,n,ind
      real*8 :: p2b,w1
      complex*16 :: cc
      real*8, allocatable :: wexa(:,:,:,:),wspl_geom(:,:,:,:,:),
     &                       wspl_asy(:,:,:,:,:),wasy(:,:,:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),
     &                           t_nu(:,:,:),aq(:,:,:,:),
     &                           q_pw(:,:,:,:,:),w_nu_nu(:,:),
     &                           q_tau(:,:,:)
	  p2b=(pi+pi)/betta_t
      i_nu=me_t*ndim3_nu+ind_nu-1
      allocate(tmp(n_pbmtm_red,n_pbmtm_red))
      allocate(tmp1(n_pbmtm_red,n_pbmtm_red))
      allocate(q_pw(n_pbmtm_red,n_pbmtm_red,2,0:n_tau/2,2))
      allocate(w_nu_nu(n_pbmtm_red,n_pbmtm_red))
      allocate(wexa(n_pbmtm_red,n_pbmtm_red,0:n_nu_exa,natom))
      allocate(wspl_geom(n_pbmtm_red,n_pbmtm_red,0:n_nu_geom+2,4,natom))
      allocate(wspl_asy(n_pbmtm_red,n_pbmtm_red,n_nu_asy+1,4,natom))
      allocate(wasy(n_pbmtm_red,n_pbmtm_red,natom))
      do iatom=1,natom
        call w_spl_loc_red(iatom,wloc(1,1,0,iatom),wexa(1,1,0,iatom),
     &                     wspl_geom(1,1,0,1,iatom),
     &                     wspl_asy(1,1,1,1,iatom),wasy(1,1,iatom))
      enddo
      nn=n_pbmt_red**2
      do ir=1,nqdiv_c
        do jatom=1,natom
          jsort=is(jatom)
          n1=n_pbmt0_red(jsort)
          jnd=iopb_red(jatom)-1
          do iatom=1,natom
            isort=is(iatom)
            n=n_pbmt0_red(isort)
            ind=iopb_red(iatom)-1
            do i0_tau=0,n_tau/2
              do il=1,2
                do it=1,2
	              q_pw(1:n,1:n1,it,i0_tau,il)=
     &              qt_pw(ind+1:ind+n,jnd+1:jnd+n1,ir,il,it,i0_tau)
                enddo
              enddo
            enddo
c -------- Transform Q(tau;nu) ---- > Q(nu';nu) ---------
            do in=1,2
              cc=(-1)**in*(0.d0,1.d0)
              allocate(aq(n,n1,0:n_tau/2,4))
              aq=(0.d0,0.d0)
              do i_tau=0,n_tau/2
                do j=1,n1
                  do i=1,n
                    aq(i,j,i_tau,1)=q_pw(i,j,1,i_tau,3-in)
     &                             +q_pw(i,j,2,i_tau,3-in)
                    aq(i,j,i_tau,2)=q_pw(i,j,1,i_tau,3-in)
     &                             -q_pw(i,j,2,i_tau,3-in)
                    aq(i,j,i_tau,2)=cc*aq(i,j,i_tau,2)
                    aq(i,j,i_tau,3)=q_pw(i,j,1,i_tau,in)
     &                             +q_pw(i,j,2,i_tau,in)
                    aq(i,j,i_tau,4)=q_pw(i,j,1,i_tau,in)
     &                             -q_pw(i,j,2,i_tau,in)
                    aq(i,j,i_tau,4)=cc*aq(i,j,i_tau,4)
                  enddo
                enddo
              enddo
              allocate(t_nu(n,n1,nomnu))
              call from_tau_nu_to_nu_nu_spl(2*n*n1,ind_nu,nomnu,t_nu,aq)
              deallocate(aq)
c ---------- Get T_q_nu and subtract static asymptotics ----------
              do i_nu1=1,nomnu
                if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
                if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)
                call interp_w_nu_re_spl(wexa(1,1,0,iatom),
     &                                  wspl_geom(1,1,0,1,iatom),
     &                                  wspl_asy(1,1,1,1,iatom),
     &                                  wasy(1,1,iatom),
     &                                  w_nu_nu,w1,n,n_pbmtm_red)
                w_nu_nu(1:n,1:n)=w_nu_nu(1:n,1:n)+vloc(1:n,1:n,iatom)
	            call zgemm('n','n',n,n1,n,(1.d0,0.d0),w_nu_nu,
     &	                   n_pbmtm_red,t_nu(1,1,i_nu1),n,(0.d0,0.d0),
     &                     tmp,n_pbmtm_red)
                if(in==1) w1=w_nu_adapt_nu(i_nu1,i_nu)
                if(in==2) w1=w_nu_adapt_nu(i_nu1,i_nu)-w_nu(i_nu)
                call interp_w_nu_re_spl(wexa(1,1,0,jatom),
     &                                  wspl_geom(1,1,0,1,jatom),
     &                                  wspl_asy(1,1,1,1,jatom),
     &                                  wasy(1,1,jatom),
     &                                  w_nu_nu,w1,n1,n_pbmtm_red)
                w_nu_nu(1:n1,1:n1)=w_nu_nu(1:n1,1:n1)
     &                            +vloc(1:n1,1:n1,jatom)
                call zgemm('n','n',n,n1,n1,(1.d0,0.d0),tmp,n_pbmtm_red,
     &                     w_nu_nu,n_pbmtm_red,(0.d0,0.d0),
     &                     t_nu(1,1,i_nu1),n)
              enddo  !! over i_nu1
c ---------------------------------------------------------------
	          if(in==1) t_nu=conjg(t_nu)
              allocate(q_tau(n,n1,0:n_tau))
              call zgemm('n','n',n*n1,n_tau+1,nomnu,(1.d0,0.d0),t_nu,
     &                   n*n1,tau_nu_from_nu_nu_1(1,0,ind_nu),
     &                   n_nu2_max,(0.d0,0.d0),q_tau,n*n1)
              deallocate(t_nu)
              if(in==1) q_tau=conjg(q_tau)
c  ---------- Add static asymptotics -----------------------
              do i0_tau=0,n_tau/2
                do it=1,2
                  i_tau=i0_tau
                  if(it==2) i_tau=n_tau-i0_tau
                  qt_pw(ind+1:ind+n,jnd+1:jnd+n1,ir,in,it,i0_tau)=
     &                   q_tau(:,:,i_tau)
                enddo
              enddo
              deallocate(q_tau)
            enddo  !! over in
          enddo   !! over iatom
        enddo  !! over jatom
      enddo  !! over ir
      deallocate(tmp,tmp1,q_pw,w_nu_nu,wexa,wspl_geom,wspl_asy,wasy)
      end
