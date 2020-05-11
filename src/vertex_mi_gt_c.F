      subroutine vertex_mi_gt_c(i_tau,ll,t_rs,grs_mi,vx,fif)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: i_tau,ll
      complex*16, intent(in) :: grs_mi(nfun_red,2*nr_full_red,0:n_tau,
     &                                 nqdiv_c),
     &                          fif(nindm_fif,maxel_red,maxel_red,
     &                              nsort),
     &                          t_rs(nr_full_red,n_pbmt_red,nqdiv_c,ll)
      complex*16, intent(inout) :: vx(nfun_red,2*nr_full_red,nqdiv_c,ll)
      integer :: iatom,isort,ind,j,il,n,i0,np,ir,il3,km,km2,ii,iii,
     &           irl,jj
      complex*16 :: cc
      real*8, allocatable :: aa(:,:,:)
      do iatom=1,natom
        isort=is(iatom)
        ind=io_lem_red(iatom)-1
        n=lfunm_red(isort)
        np=n_pbmt0_red(isort)
	    i0=iopb_red(iatom)-1
	    allocate(aa(n,n,nr_full_red))
        do ir=1,nqdiv_c
          do il=1,ll
            il3=3-il
            if(ll==1) il3=1
            aa=(0.d0,0.d0)
            do km=1,n
              do km2=1,n
                do ii=1,nind_fif(km2,km,isort)
                  iii=i0+ind_fif(ii,km2,km,isort)
                  cc=conjg(fif(ii,km2,km,isort))
                  aa(km2,km,:)=aa(km2,km,:)+cc*t_rs(:,iii,ir,il3)
                enddo
              enddo
            enddo
            do j=1,nr_full_red
              do irl=1,2
                jj=(irl-1)*nr_full_red+j
                do km=1,n
                  cc=(0.d0,0.d0)
                  do km2=1,n
                    cc=cc+grs_mi(ind+km2,jj,i_tau,ir)*aa(km2,km,j)
                  enddo
                  if(ll==2) cc=cc*eqr(ir)
                  vx(ind+km,jj,ir,il)=vx(ind+km,jj,ir,il)+cc
                enddo
              enddo
            enddo   !! over j
          enddo    !! over il
        enddo   !! over ir
        deallocate(aa)
      enddo  !! over iatom
      end
