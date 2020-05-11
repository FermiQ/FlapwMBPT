      subroutine map_pbr_scalapack
      use manager_mod
      use parallel_mod
      use scalapack_mod
      implicit none
      integer :: m,i,j
      nrow_b=int(sqrt(real(nproc_b)))
1     ncol_b=nproc_b/nrow_b
      if(ncol_b*nrow_b/=nproc_b) then
        nrow_b=nrow_b+1
        goto 1
      endif
      if(maswrk) write(iun,*)' MAP_PBR_SCALAPACK: ',' nrow_b = ',nrow_b,
     &                       ' ncol_b = ',ncol_b
      allocate(map_pbr(nrow_b,ncol_b))
      map_pbr=0
      do m=0,nproc_t*nproc_k-1
        if(me_t_k==m) then
          i=mod(me_b,nrow_b)+1
          j=me_b/nrow_b+1
          map_pbr(i,j)=me
        endif
      enddo
      if(nproc_b/=1) call igop(map_pbr,nrow_b*ncol_b,'  +',comm_b)
      lsclp=.false.
#ifdef SCLP
      if(nproc_b>1) lsclp=.true.
#endif
      end
