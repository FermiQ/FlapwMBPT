      subroutine output_vrt_wi
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,ispin,nk,i,ka,i_len,info,j
      complex*16 :: cc
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: tmv(:,:),vuni(:,:)
      cc=dcmplx(0.d0,-(pi+pi)/betta_t)
      if(allocated(vrt_scale)) deallocate(vrt_scale)
      allocate(vrt_scale(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_VRT_WI_0')
        open(4,file=allfile(1:i_len)//'_VRT_WI_1')
      endif
      do ispin=1,nspin
        do k=1,npnt_c
          ka=k_a_from_c(k)
          nk=n_low_bnd(ka,ispin)
          allocate(tmv(nk,nk))
          do j=1,nk
            do i=1,nk
              tmv(i,j)=(sig_c_wi(i,j,k,1,ispin)-sig_c_wi(i,j,k,0,ispin))
     &                /cc
            enddo
          enddo
          if(maswrk) then
            write(3,*)' Ispin K = ',ispin,k
            write(4,*)' Ispin K = ',ispin,k
            do i=1,nk
              write(3,'(i4,2(1x,f11.6))')i,real(vrt_wi(i,i,k,0,ispin)),
     &                                   real(z_wi(i,i,k,ispin))
              write(4,'(i4,4(1x,f11.6))')i,vrt_wi(i,i,k,1,ispin),
     &                                   tmv(i,i)
            enddo
          endif
          allocate(ipiv(nk))
          allocate(vuni(nk,nk))
          if(metal) then
            vuni=vrt_wi(1:nk,1:nk,k,1,ispin)
            vrt_scale(1:nk,1:nk,k,ispin)=tmv
          else
            vuni=0.5d0*(vrt_wi(1:nk,1:nk,k,0,ispin)
     &                 +conjg(transpose(vrt_wi(1:nk,1:nk,k,0,ispin))))
            vrt_scale(1:nk,1:nk,k,ispin)=z_wi(1:nk,1:nk,k,ispin)
          endif
          deallocate(tmv)
          call zgesv(nk,nk,vuni,nk,ipiv,vrt_scale(1,1,k,ispin),
     &               nbndf_bnd,info)
          deallocate(ipiv,vuni)
        enddo   !! over k
      enddo   !! over ispin
      if(maswrk) then
        close(3)
        close(4)
      endif
      end
