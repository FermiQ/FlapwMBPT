      subroutine t_matr_sf_gen(tmatr,g_rr_tau,wst4loc,ch,nt)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: nt
      real*8, intent(in) :: g_rr_tau(nfun_red,nfun_red,nqdiv_c,2,
     &                               ndim3_tau,nspin),
     &                      wst4loc(maxel_red**2,maxel_red**2,natom,
     &                              nspin,nspin)
      complex*16, intent(out) :: tmatr(n_lem2_red,n_lem2_red,ndim3_nu,
     &                                 ndimc_kk(me_k+1),nspin,nspin)
      integer :: mx2,iatom,isort,ind2,n,nn,ispin,l,i,ind_nu,ind_q,jspin,
     &           k,kl,ij,j,ii,jj
      complex*16, allocatable :: vc(:,:,:,:,:),a_ph(:,:,:,:,:,:)
      mx2=maxel_red**2
c --- Reordering W according to PP or PH type -------------------------
      allocate(vc(mx2,mx2,natom,nspin,nspin))
      do jspin=1,nspin
        do ispin=1,nspin
          do iatom=1,natom
            isort=is(iatom)
            n=lfunm_red(isort)
            kl=0
            do l=1,n
              do k=1,n
                kl=kl+1
                ij=0
                do j=1,n
                  do i=1,n
                    ij=ij+1
                    if(ch=='pp') then
                      ii=(k-1)*n+i
                      jj=(j-1)*n+l
                      vc(ij,kl,iatom,ispin,jspin)=
     &                    wst4loc(ii,jj,iatom,ispin,jspin)
                    else if(ch=='ph') then
                      ii=(k-1)*n+i
                      jj=(l-1)*n+j
                      vc(ij,kl,iatom,ispin,jspin)=
     &                    wst4loc(ii,jj,iatom,ispin,jspin)
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
c --- 4-point particle-particle Response function ------------------
      do jspin=1,nspin
        do ispin=1,nspin
          call l_q_nu_g2_w1(ispin,jspin,tmatr(1,1,1,1,ispin,jspin),
     &                      g_rr_tau,ch)
        enddo
      enddo
c --- Auxiliary matrices A=WL ----------------------
      allocate(a_ph(n_lem2_red,n_lem2_red,ndim3_nu,
     &              ndimc_kk(me_k+1),nspin,nspin))
      do jspin=1,nspin
        do ispin=1,nspin
          do ind_nu=1,ndim3_nu
            do ind_q=1,ndimc_kk(me_k+1)
              do iatom=1,natom
                isort=is(iatom)
                n=lfunm_red(isort)
                nn=n*n
                ind2=io_lem2_red(iatom)
                call zgemm('n','n',nn,n_lem2_red,nn,(-1.d0,0.d0),
     &                     vc(1,1,iatom,ispin,jspin),mx2,
     &                     tmatr(ind2,1,ind_nu,ind_q,ispin,jspin),
     &                     n_lem2_red,(0.d0,0.d0),
     &                     a_ph(ind2,1,ind_nu,ind_q,ispin,jspin),
     &                     n_lem2_red)
              enddo
            enddo
          enddo
        enddo
      enddo
c ---   T-matrix in (q;nu) variables -------------------------------
      do jspin=1,nspin
        do ispin=1,nspin
          call t_matr_inf(vc(1,1,1,ispin,jspin),
     &                    a_ph(1,1,1,1,ispin,jspin),nt,
     &                    tmatr(1,1,1,1,ispin,jspin))
        enddo
      enddo
      deallocate(a_ph,vc)
      call timel('* T_MATR_SF_GEN finished ***********')
      end
