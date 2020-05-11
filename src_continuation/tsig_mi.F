      subroutine tsig_ii(g_rs_ii,k_rs,t_rs,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      real*8, intent(in) :: g_rs_ii(nr_red_red,nr_full_red,0:n_tau,
     &                              nqdiv_c)
      complex*16, intent(in) :: k_rs(nr_full_red,nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(nr_full_red,nr_full_red,nqdiv_c)
      integer :: i,j,j0,i0,ig,ir,ir0
      real*8 :: v(3),v1(3)
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
            call zone1_number(v1,gbas,ndiv_c,ir0)
            j0=ip_r_red(j,ig)
            t_rs(i,j,ir)=-g_rs_ii(i0,j0,i1_tau,ir0)*k_rs(i,j,ir)
          enddo   !! over j
        enddo   !! over ir
      enddo   !! over i
      end
      
      
      
      subroutine tsig_ii_c(g_rs_ii,k_rs,t_rs,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      complex*16, intent(in) :: g_rs_ii(2*nr_red_red,2*nr_full_red,
     &                                  0:n_tau,nqdiv_c)
      complex*16, intent(in) :: k_rs(2*nr_full_red,nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(nr_full_red,2*nr_full_red,nqdiv_c)
      integer :: i,j,j0,i0,ig,ir,ir0,irl,jrl,i1,j1,is0,js0,is1,js,ir2
      real*8 :: v(3),v1(3)
      complex*16 :: cc,gc,cc1
      t_rs=(0.d0,0.d0)
      do i=1,nr_full_red
        i0=ir_ref_red(i)
        ig=ir_group_red(i)
        do ir=1,nqdiv_c
          do j=1,nr_full_red
            v=rr_coarse(:,ir)+tshift_r_red(:,i,ig)-tshift_r_red(:,j,ig)
            call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
            call zone1_number(v1,gbas,ndiv_c,ir0)
            j0=ip_r_red(j,ig)
            do jrl=1,2
              j1=(jrl-1)*nr_full_red+j
              do js=1,2
                js0=(js-1)*nr_full_red+j0
                cc1=conjg(uj((js-1)*2+jrl,ig))
                do irl=1,2
                  i1=(irl-1)*nr_full_red+i
                  ir2=irl-2
                  do is1=1,2
                    is0=(is1-1)*nr_red_red+i0
                    cc=uj(is1+is1+ir2,ig)*cc1
                    gc=conjg(g_rs_ii(is0,js0,i1_tau,ir0))
                    t_rs(i,j1,ir)=t_rs(i,j1,ir)-cc*gc*k_rs(i1,j,ir)
                  enddo
                enddo
              enddo
            enddo
          enddo   !! over j
        enddo   !! over ir
      enddo   !! over i
      end
      
      
      subroutine tsig_mi(g_rs_mi,k_rs,t_rs,fif,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      real*8, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,nsort),
     &                      g_rs_mi(nfun_red,nr_full_red,0:n_tau,
     &                              nqdiv_c)
      complex*16, intent(in) :: k_rs(nfun_red,nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(n_pbmt_red,nr_full_red,nqdiv_c)
      integer :: mx2,iatom,isort,ind,i,j,ii,nn,i0,n0,ij,n00,i1,ir,iii
      real*8 :: tt
      complex*16, allocatable :: dd(:,:),ft(:,:)
      mx2=maxel_red**2
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n0=lfunm_red(isort)
          n00=n_pbmt0_red(isort)
	      i0=iopb_red(iatom)-1
          nn=n0*n0
	      allocate(ft(nr_full_red,nn))
          ij=0
	      do j=1,n0
	        do i=1,n0
              ij=ij+1
              do ii=1,nr_full_red
                ft(ii,ij)=g_rs_mi(ind+i,ii,i1_tau,ir)*k_rs(ind+j,ii,ir)
              enddo
            enddo
          enddo
	      allocate(dd(nr_full_red,n00))
          dd=(0.d0,0.d0)
          ij=0
          do j=1,n0
            do i=1,n0
              ij=ij+1
              do ii=1,nind_fif(i,j,isort)
                iii=ind_fif(ii,i,j,isort)
                tt=fif(ii,i,j,isort)
                dd(:,iii)=dd(:,iii)+tt*ft(:,ij)
              enddo
            enddo
          enddo
          do j=1,nr_full_red
            do i=1,n00
              i1=i0+i
              t_rs(i1,j,ir)=-dd(j,i)
            enddo
          enddo
          deallocate(ft,dd)
        enddo  !! over iatom
      enddo   !! over ir
      end
      
      
      subroutine tsig_mi_c(g_rs_mi,k_rs,t_rs,fif,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      complex*16, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,
     &                              nsort),
     &                          g_rs_mi(nfun_red,2*nr_full_red,0:n_tau,
     &                                  nqdiv_c),
     &                          k_rs(nfun_red,nr_full_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(n_pbmt_red,2*nr_full_red,nqdiv_c)
      integer :: mx2,iatom,isort,ind,i,j,ii,nn,i0,n0,ij,n00,i1,ir,
     &           irl,iii
      complex*16, allocatable :: dd(:,:),ft(:,:)
      mx2=maxel_red**2
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n0=lfunm_red(isort)
          n00=n_pbmt0_red(isort)
	      i0=iopb_red(iatom)-1
          nn=n0*n0
	      allocate(ft(nn,2*nr_full_red))
          ij=0
	      do j=1,n0
	        do i=1,n0
              ij=ij+1
              do ii=1,nr_full_red
                do irl=1,2
                  iii=(irl-1)*nr_full_red+ii
	              ft(ij,iii)=conjg(g_rs_mi(ind+i,iii,i1_tau,ir))
     &                      *k_rs(ind+j,ii,ir)
                enddo
              enddo
            enddo
          enddo
	      allocate(dd(n00,2*nr_full_red))
          dd=(0.d0,0.d0)
          ij=0
          do j=1,n0
            do i=1,n0
              ij=ij+1
              do ii=1,nind_fif(i,j,isort)
                iii=ind_fif(ii,i,j,isort)
                dd(iii,:)=dd(iii,:)+fif(ii,i,j,isort)*ft(ij,:)
              enddo
            enddo
          enddo
          do j=1,2*nr_full_red
            do i=1,n00
              i1=i0+i
              t_rs(i1,j,ir)=-dd(i,j)
            enddo
          enddo
          deallocate(ft,dd)
        enddo  !! over iatom
      enddo   !! over ir
      end
      
      
      
      
      subroutine tsig_im(g_rs_im,k_rs,t_rs,fif,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      real*8, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,nsort),
     &                      g_rs_im(nr_full_red,nfun_red,0:n_tau,
     &                              nqdiv_c)
      complex*16, intent(in) :: k_rs(nr_full_red,n_pbmt_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(nr_full_red,nfun_red,nqdiv_c)
      integer :: iatom,isort,ind,i,ii,indp,n,np,ir,km1,km2,
     &           iii,jj,jj0
      real*8, allocatable :: aa(:,:)
      t_rs=(0.d0,0.d0)
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
	      indp=iopb_red(iatom)-1
          allocate(aa(nr_full_red,np))
          do km1=1,n
            aa=0.d0
            do km2=1,n
              do ii=1,nind_fif(km1,km2,isort)
                iii=ind_fif(ii,km1,km2,isort)
                aa(:,iii)=aa(:,iii)+fif(ii,km1,km2,isort)
     &                             *g_rs_im(:,ind+km2,i1_tau,ir)
              enddo
            enddo
            ii=ind+km1
            do jj=1,np
              jj0=indp+jj
              do i=1,nr_full_red
                t_rs(i,ii,ir)=t_rs(i,ii,ir)-k_rs(i,jj0,ir)*aa(i,jj)
              enddo
            enddo
          enddo
          deallocate(aa)
        enddo  !! over iatom
      enddo   !! over ir
      end
      
      
      
      
      subroutine tsig_im_c(g_rs_im,k_rs,t_rs,fif,i1_tau)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      complex*16, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,
     &                              nsort),
     &                          g_rs_im(2*nr_full_red,nfun_red,0:n_tau,
     &                                  nqdiv_c),
     &                          k_rs(2*nr_full_red,n_pbmt_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(nr_full_red,nfun_red,nqdiv_c)
      integer :: iatom,isort,ind,i,ii,indp,n,np,ir,km1,km2,
     &           iii,jj,jj0,irl,ir0,m1,m2
      real*8, allocatable :: aa(:,:)
      t_rs=(0.d0,0.d0)
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
	      indp=iopb_red(iatom)-1
          allocate(aa(nr_full_red,np))
          do km1=1,n
            do irl=1,2
              ir0=(irl-1)*nr_full_red
              m1=ir0+1
              m2=ir0+nr_full_red
              aa=0.d0
              do km2=1,n
                do ii=1,nind_fif(km1,km2,isort)
                  iii=ind_fif(ii,km1,km2,isort)
                  aa(:,iii)=aa(:,iii)+fif(ii,km1,km2,isort)
     &                      *conjg(g_rs_im(m1:m2,ind+km2,i1_tau,ir))
                enddo
              enddo
              ii=ind+km1
              do jj=1,np
                jj0=indp+jj
                do i=1,nr_full_red
                  t_rs(i,ii,ir)=t_rs(i,ii,ir)-k_rs(ir0+i,jj0,ir)
     &                                       *aa(i,jj)
                enddo
              enddo
            enddo
          enddo
          deallocate(aa)
        enddo  !! over iatom
      enddo   !! over ir
      end
