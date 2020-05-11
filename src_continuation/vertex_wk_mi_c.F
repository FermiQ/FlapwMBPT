      subroutine vertex_wk_mi_c(wmi,k_rs,vx,fif)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: wmi(n_pbmt_red,nr_full_red,nqdiv_c)
      complex*16, intent(in) :: k_rs(nfun_red,2*nr_full_red,nqdiv_c),
     &                          fif(nindm_fif,maxel_red,maxel_red,nsort)
      complex*16, intent(inout) :: vx(nfun_red,2*nr_full_red,nqdiv_c)
      integer :: iatom,isort,ind,i,j,n,i0,ir,irl,ii,iii,irr
      complex*16 :: tt
      complex*16, allocatable :: bb(:,:),k0(:,:),aa(:,:,:),a0(:,:)
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          i0=iopb_red(iatom)-1
          allocate(aa(n,n,nr_full_red))
          aa=(0.d0,0.d0)
          do j=1,n
            do i=1,n
              do ii=1,nind_fif(i,j,isort)
                iii=i0+ind_fif(ii,i,j,isort)
                tt=fif(ii,i,j,isort)
                aa(i,j,:)=aa(i,j,:)+wmi(iii,:,ir)*tt
              enddo
            enddo
          enddo
          allocate(k0(2,n))
          allocate(a0(n,n))
          allocate(bb(2,n))
          do irr=1,nr_full_red
            do irl=1,2
              ii=(irl-1)*nr_full_red
              do i=1,n
                k0(irl,i)=k_rs(ind+i,ii+irr,ir)
              enddo
            enddo
            a0=transpose(aa(:,:,irr))
            call zgemm('n','n',2,n,n,(1.d0,0.d0),k0,2,a0,n,(0.d0,0.d0),
     &                 bb,2)
            do irl=1,2
              ii=(irl-1)*nr_full_red+irr
              do i=1,n
                vx(ind+i,ii,ir)=vx(ind+i,ii,ir)+bb(irl,i)
              enddo
            enddo
          enddo   !! over irr
          deallocate(aa,k0,a0,bb)
        enddo  !! over iatom
      enddo   !! over ir
      end
