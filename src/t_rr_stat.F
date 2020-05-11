      subroutine t_rr_stat(qt_pw,vloc,ind_nu,nomnu)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ind_nu,nomnu
      real*8, intent(in) :: vloc(n_pbmtm_red,n_pbmtm_red,natom)
      complex*16, intent(inout) :: qt_pw(n_pbmt_red,n_pbmt_red,
     &                                   nqdiv_c,2,2,0:n_tau/2)
      integer :: j,i,nn,i_nu1,in,i_tau,i_nu,i0_tau,it,il,iatom,isort,
     &           jatom,jsort,ir,n1,jnd,n,ind
      real*8 :: p2b
      complex*16 :: cc
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),
     &                           t_nu(:,:,:),aq(:,:,:,:),
     &                           q_pw(:,:,:,:,:),q_tau(:,:,:)
      p2b=(pi+pi)/betta_t
      i_nu=me_t*ndim3_nu+ind_nu-1
      allocate(tmp(n_pbmtm_red,n_pbmtm_red))
      allocate(tmp1(n_pbmtm_red,n_pbmtm_red))
      allocate(q_pw(n_pbmtm_red,n_pbmtm_red,2,0:n_tau/2,2))
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
                call zgemm('n','n',n,n1,n,(1.d0,0.d0),vloc(1,1,iatom),
     &                     n_pbmtm_red,t_nu(1,1,i_nu1),n,(0.d0,0.d0),
     &                     tmp,n_pbmtm_red)
                call zgemm('n','n',n,n1,n1,(1.d0,0.d0),tmp,n_pbmtm_red,
     &                     vloc(1,1,jatom),n_pbmtm_red,(0.d0,0.d0),
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
      deallocate(tmp,tmp1,q_pw)
      end
