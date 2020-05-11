      subroutine vertex_wk_loc_r(wpb,k_rs,vx,fif,isp,jsp)
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
     &                      wpb(n_pbmtm_red,n_pbmtm_red,natom)
      complex*16, intent(in) :: k_rs(maxel_red,maxel_red,natom)
      complex*16, intent(inout) :: vx(maxel_red,maxel_red,natom)
      integer :: iatom,isort,i,j,ii,iii,n,np
      real*8 :: tt
      real*8, allocatable :: aa(:,:,:),a0(:,:,:)
      complex*16, allocatable :: bb(:,:,:),dd(:,:),kr0(:,:)
      do iatom=1,natom
        isort=is(iatom)
        n=lfunm_red(isort)
        np=n_pbmt0_red(isort)
        allocate(aa(np,n,n))
        aa=0.d0
        do j=1,n
          do i=1,n
            do ii=1,nind_fif(i,j,isort)
              iii=ind_fif(ii,i,j,isort)
              tt=fif(ii,i,j,isort,isp)
              aa(:,i,j)=aa(:,i,j)+wpb(iii,1:np,iatom)*tt
            enddo
          enddo
        enddo
        allocate(kr0(n,n))
        do i=1,n
          do j=1,n
            kr0(j,i)=k_rs(i,j,iatom)
          enddo
        enddo
        allocate(a0(n,np,n))
	    do j=1,n
          do ii=1,np
	        do i=1,n
	          a0(i,ii,j)=aa(ii,j,i)
            enddo
	      enddo
	    enddo
	    allocate(bb(n,np,n))
	    call dgemm('n','n',2*n,np*n,n,1.d0,kr0,2*n,a0,n,0.d0,bb,2*n)
	    allocate(dd(n,n))
        dd=(0.d0,0.d0)
        do j=1,n
          do i=1,n
            do ii=1,nind_fif(i,j,isort)
              iii=ind_fif(ii,i,j,isort)
              tt=fif(ii,i,j,isort,jsp)
              dd(:,j)=dd(:,j)+bb(i,iii,:)*tt
            enddo
          enddo
        enddo
        do j=1,n
          do i=1,n
            vx(i,j,iatom)=dd(i,j)
          enddo
        enddo
        deallocate(kr0,a0,bb,dd)
        deallocate(aa)
      enddo   !! over iatom
      end
