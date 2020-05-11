#ifdef SCLP
      subroutine sclp_pw_prep
      use atom_mod
      use scalapack_mod
      use solid_mod
      implicit none
      integer :: myrow,mycol,k,numroc,n,i,ii,info
      nblock_pw=64
      call blacs_get(0,0,sctxt_pw)
      call blacs_gridmap(sctxt_pw,map_pbr,nrow_b,nrow_b,ncol_b)
      call blacs_gridinfo(sctxt_pw,nrow_b,ncol_b,myrow,mycol)
      allocate(ncol_pw(npnt))
      allocate(nrow_pw(npnt))
      do k=1,npnt
        n=n_pbmt+nplwgw(k)
        nrow_pw(k)=numroc(n,nblock_pw,myrow,0,nrow_b)
        ncol_pw(k)=numroc(n,nblock_pw,mycol,0,ncol_b)
      enddo
      ncol_pw_max=maxval(ncol_pw)
      nrow_pw_max=maxval(nrow_pw)
      allocate(refrow_pw(nrow_pw_max,npnt))
      allocate(refcol_pw(ncol_pw_max,npnt))
      allocate(desc_pw(9,npnt))
      do k=1,npnt
        n=n_pbmt+nplwgw(k)
        do i=1,ncol_pw(k)
          ii=(i-1)/nblock_pw
          ii=(ii*ncol_b+mycol)*nblock_pw+i-ii*nblock_pw
          refcol_pw(i,k)=ii
        enddo
        do i=1,nrow_pw(k)
          ii=(i-1)/nblock_pw
          ii=(ii*nrow_b+myrow)*nblock_pw+i-ii*nblock_pw
          refrow_pw(i,k)=ii
        enddo
        call descinit(desc_pw(1,k),n,n,nblock_pw,nblock_pw,0,0,sctxt_pw,
     &                nrow_pw(k),info)
      enddo
      end
#endif
