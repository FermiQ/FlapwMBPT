      subroutine vertex_wk_loc_c(wpb,k_rs,vx,fif)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      real*8, intent(in) :: wpb(n_pbmtm_red,n_pbmtm_red,natom)
      complex*16, intent(in) :: k_rs(maxel_red,maxel_red,natom),
     &                          fif(nindm_fif,maxel_red,maxel_red,nsort)
      complex*16, intent(inout) :: vx(maxel_red,maxel_red,natom)
      integer :: iatom,isort,i,j,ii,iii,n,np
      complex*16 :: tt
      complex*16, allocatable :: bb(:,:,:),dd(:,:),kr0(:,:),aa(:,:,:),
     &                           a0(:,:,:)
      do iatom=1,natom
        isort=is(iatom)
        n=lfunm_red(isort)
        np=n_pbmt0_red(isort)
        allocate(aa(np,n,n))
        aa=(0.d0,0.d0)
        do j=1,n
          do i=1,n
            do ii=1,nind_fif(i,j,isort)
              iii=ind_fif(ii,i,j,isort)
              tt=fif(ii,i,j,isort)
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
	  call zgemm('n','n',n,np*n,n,(1.d0,0.d0),kr0,n,a0,n,(0.d0,0.d0),
     &             bb,n)
        allocate(dd(n,n))
        dd=(0.d0,0.d0)
        do j=1,n
          do i=1,n
            do ii=1,nind_fif(i,j,isort)
              iii=ind_fif(ii,i,j,isort)
              tt=fif(ii,i,j,isort)
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
