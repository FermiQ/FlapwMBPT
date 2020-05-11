      subroutine output_vrt_uni(lambda_stat,lambda_dyn)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: lambda_dyn(nbndf_bnd,nbndf_bnd,2,
     &                                 ndim3_omega,npnt_c,nspin)
      complex*16, intent(in) :: lambda_stat(nbndf_bnd,nbndf_bnd,npnt_c,
     &                                      nspin)
      logical :: scl
      integer :: i_omega,k,ispin,nk,i,ka,ind_omega,nb,i_len,info
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: pr(:,:),tmv(:,:),vuni(:,:)
      scl=.false.
      if(scale_psi.and.iter==iter_dft+iter_hf+iter_gw+iter_qp+1)
     &  scl=.true.
      if(scl) then
        if(allocated(vrt_scale)) deallocate(vrt_scale)
        allocate(vrt_scale(nbndf_bnd,nbndf_bnd,npnt_c,nspin))
      endif
      nb=min(6,minval(n_low_bnd))
      allocate(pr(0:n_omega,nb))
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_Im_'//'_VRT_UNI')
        open(4,file=allfile(1:i_len)//'_Re_'//'_VRT_UNI')
      endif
      if(scl) vrt_scale=z_wi
      do ispin=1,nspin
        do k=1,npnt_c
          ka=k_a_from_c(k)
          nk=n_low_bnd(ka,ispin)
          allocate(tmv(nk,nk))
          pr=(0.d0,0.d0)
          if(scl) then
            allocate(vuni(nk,nk))
            vuni=(0.d0,0.d0)
          endif
          do ind_omega=1,ndim3_omega
            i_omega=me_t*ndim3_omega+ind_omega-1
            call ferm_unpack_omega(tmv,
     &                             lambda_dyn(1,1,1,ind_omega,k,ispin),
     &                             nk,nbndf_bnd,nk)
            do i=1,nb
              pr(i_omega,i)=tmv(i,i)+lambda_stat(i,i,k,ispin)
            enddo
            if(i_omega==0) then
              tmv=tmv+lambda_stat(1:nk,1:nk,k,ispin)
              if(scl) vuni=0.5d0*(tmv+transpose(conjg(tmv)))
            endif
          enddo  !! over i_omega
          deallocate(tmv)
          if(nproc_t/=1) then
            call dgop(pr,2*(n_omega+1)*nb,'  +',comm_t)
            if(scl) call dgop(vuni,2*nk*nk,'  +',comm_t)
          endif
          pr=pr*sqrt(amega)
          if(scl) vuni=vuni*sqrt(amega)
          if(maswrk) then
            write(3,*)' Ispin K = ',ispin,k
            write(4,*)' Ispin K = ',ispin,k
            do i_omega=0,n_omega
              write(3,'(f11.4,6(1x,f11.6))')w_omega(i_omega),
     &                     (imag(pr(i_omega,i)),i=1,nb)
              write(4,'(f11.4,6(1x,f11.6))')w_omega(i_omega),
     &                     (real(pr(i_omega,i)),i=1,nb)
            enddo
          endif
          if(scl) then
            allocate(ipiv(nk))
            call zgesv(nk,nk,vuni,nk,ipiv,vrt_scale(1,1,k,ispin),
     &                 nbndf_bnd,info)
            deallocate(ipiv,vuni)
          endif
        enddo   !! over k
      enddo   !! over ispin
      deallocate(pr)
      if(maswrk) then
        close(3)
        close(4)
      endif
      end
