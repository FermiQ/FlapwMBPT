      subroutine vertex_wk_mm_r(wpb,k_rs,vx,fif,isp,jsp)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: isp,jsp
      real*8, intent(in) :: fif(nindm_fif,maxel_red,maxel_red,nsort,
     &                          nspin),
     &                      wpb(n_pbmt_red,n_pbmt_red,nqdiv_c)
      complex*16, intent(in) :: k_rs(nfun_red,nfun_red,nqdiv_c)
      complex*16, intent(inout) :: vx(nfun_red,nfun_red,nqdiv_c)
      integer :: jatom,jsort,jnd,iatom,isort,ind,i,j,j0,n1,i0,
     &           j1,i1,ir,ii,iii,n,np,np1
      real*8 :: tt
      real*8, allocatable :: aa(:,:,:),a0(:,:,:)
      complex*16, allocatable :: bb(:,:,:),dd(:,:),kr0(:,:)
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
	      i0=iopb_red(iatom)-1
          allocate(aa(n_pbmt_red,n,n))
          aa=0.d0
          do j=1,n
            do i=1,n
              do ii=1,nind_fif(i,j,isort)
                iii=i0+ind_fif(ii,i,j,isort)
                tt=fif(ii,i,j,isort,isp)
                aa(:,i,j)=aa(:,i,j)+wpb(iii,:,ir)*tt
              enddo
            enddo
          enddo
          do jatom=1,natom
            jsort=is(jatom)
            jnd=io_lem_red(jatom)-1
            n1=lfunm_red(jsort)
            np1=n_pbmt0_red(jsort)
	        j0=iopb_red(jatom)-1
	        allocate(kr0(n1,n))
            do i=1,n
              do j=1,n1
                kr0(j,i)=k_rs(ind+i,jnd+j,ir)
              enddo
            enddo
            allocate(a0(n,np1,n))
	        do j=1,n
              do ii=1,np1
	            do i=1,n
	              a0(i,ii,j)=aa(j0+ii,j,i)
                enddo
	          enddo
	        enddo
	        allocate(bb(n1,np1,n))
	        call dgemm('n','n',2*n1,np1*n,n,1.d0,kr0,2*n1,a0,n,0.d0,
     &                 bb,2*n1)
	        allocate(dd(n,n1))
            dd=(0.d0,0.d0)
            do j=1,n1
              do i=1,n1
                do ii=1,nind_fif(i,j,jsort)
                  iii=ind_fif(ii,i,j,jsort)
                  tt=fif(ii,i,j,jsort,jsp)
                  dd(:,j)=dd(:,j)+bb(i,iii,:)*tt
                enddo
              enddo
            enddo
            do j=1,n1
              j1=jnd+j
              do i=1,n
                i1=ind+i
                vx(i1,j1,ir)=dd(i,j)
              enddo
            enddo
            deallocate(kr0,a0,bb,dd)
          enddo  !! over jatom
          deallocate(aa)
        enddo   !! over iatom
      enddo   !! over ir
      end
