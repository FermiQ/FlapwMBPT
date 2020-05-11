      subroutine tsig_mm(g_rs_mm,k_rs,t_rs,fif,i1_tau)
  	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i1_tau
      real*8, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,nsort),
     &                      g_rs_mm(nfun_red,nfun_red,0:n_tau,nqdiv_c)
      complex*16, intent(in) :: k_rs(nfun_red,n_pbmt_red,nqdiv_c)
      complex*16, intent(out) :: t_rs(n_pbmt_red,nfun_red,nqdiv_c)
      integer :: mx2,jatom,jsort,jnd,iatom,isort,ind,i,j,nn,j0,n1,i0,
     &           n0,n11,n00,j1,i1,ir,ii,iii
      real*8 :: tt
      real*8, allocatable :: gg(:,:)
      complex*16, allocatable :: aa(:,:,:),bb(:,:,:),dd(:,:),ft(:,:,:)
      mx2=maxel_red**2
      do ir=1,nqdiv_c
        do jatom=1,natom
          jsort=is(jatom)
          jnd=io_lem_red(jatom)-1
          n1=lfunm_red(jsort)
          n11=n_pbmt0_red(jsort)
	      j0=iopb_red(jatom)-1
          allocate(aa(nfun_red,n1,n1))
          aa=(0.d0,0.d0)
          do j=1,n1
            do i=1,n1
              do ii=1,nind_fif(i,j,jsort)
                iii=j0+ind_fif(ii,i,j,jsort)
                tt=fif(ii,i,j,jsort)
                aa(:,i,j)=aa(:,i,j)+k_rs(:,iii,ir)*tt
              enddo
            enddo
          enddo
          do iatom=1,natom
            isort=is(iatom)
            ind=io_lem_red(iatom)-1
            n0=lfunm_red(isort)
            n00=n_pbmt0_red(isort)
	        i0=iopb_red(iatom)-1
            nn=n0*n1
	        allocate(ft(n1,n0,n1))
	        do j=1,n1
	          do i=1,n1
	            ft(i,:,j)=aa(ind+1:ind+n0,i,j)
	          enddo
	        enddo
	        allocate(bb(n1,n0,n0))
            allocate(gg(n0,n1))
            gg=g_rs_mm(ind+1:ind+n0,jnd+1:jnd+n1,i1_tau,ir)
	        call dgemm('n','t',2*nn,n0,n1,1.d0,ft,2*nn,gg,n0,0.d0,
     &                 bb,2*nn)
	        allocate(dd(n1,n00))
            dd=(0.d0,0.d0)
            do j=1,n0
              do i=1,n0
                do ii=1,nind_fif(i,j,isort)
                  iii=ind_fif(ii,i,j,isort)
                  tt=fif(ii,i,j,isort)
                  dd(:,iii)=dd(:,iii)+bb(:,j,i)*tt
                enddo
              enddo
            enddo
            do j=1,n1
              j1=jnd+j
              do i=1,n00
                i1=i0+i
                t_rs(i1,j1,ir)=-dd(j,i)
              enddo
            enddo
            deallocate(ft,bb,dd,gg)
          enddo  !! over iatom
          deallocate(aa)
        enddo   !! over jatom
      enddo   !! over ir
      end
