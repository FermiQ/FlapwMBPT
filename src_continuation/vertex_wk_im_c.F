      subroutine vertex_wk_im_c(wim,k_rs,vx,fif)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: wim(nr_full_red,n_pbmt_red,nqdiv_c)
      complex*16, intent(in) :: k_rs(2*nr_full_red,nfun_red,nqdiv_c),
     &                          fif(nindm_fif,maxel_red,maxel_red,nsort)
      complex*16, intent(inout) :: vx(2*nr_full_red,nfun_red,nqdiv_c)
      integer :: iatom,isort,ind,i,j,n,np,i0,ir,irl,jj,i1,i2,i3,ii,iii
      complex*16 :: tt
      complex*16, allocatable :: a(:)
      allocate(a(nr_full_red))
      do ir=1,nqdiv_c
        do iatom=1,natom
          isort=is(iatom)
          ind=io_lem_red(iatom)-1
          n=lfunm_red(isort)
          np=n_pbmt0_red(isort)
          i0=iopb_red(iatom)-1
          do j=1,n
            jj=ind+j
            do i=1,n
              a=(0.d0,0.d0)
              do ii=1,nind_fif(i,j,isort)
                iii=i0+ind_fif(ii,i,j,isort)
                tt=fif(ii,i,j,isort)
                a(:)=a(:)+tt*wim(:,iii,ir)
              enddo
              ii=ind+i
              do irl=1,2
                i1=(irl-1)*nr_full_red
                do i2=1,nr_full_red
                  i3=i1+i2
                  vx(i3,jj,ir)=vx(i3,jj,ir)+a(i2)*k_rs(i3,ii,ir)
                enddo
              enddo
            enddo
          enddo   !! over j
        enddo  !! over iatom
      enddo   !! over ir
      deallocate(a)
      end
