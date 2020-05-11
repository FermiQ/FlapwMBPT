      subroutine wsq0_ins_r(p00,p_head,p_wing,scal,ncol,nrow,refcol,
     &                      refrow,nblock,desc,nne,w_s,n1,n2,eout,e1out,
     &                      ind_nu,i_nu,vsq)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
#ifdef SCLP
      integer :: i1,j1
#endif
      logical, intent(in) :: scal
      integer, intent(in) :: nne,ncol,nrow,refcol(ncol),refrow(nrow),
     &                       nblock,desc(9),n1,n2,ind_nu,i_nu
      real*8, intent(in) :: vsq(nne)
      real*8, intent(in) :: p_head(3,3,ndim3_nu),
     &                      p_wing(3,nne,ndim3_nu)
      real*8, intent(out) :: p00
      real*8, intent(inout) :: w_s(n1,n2)
      complex*16, intent(out) :: eout(3,3,0:n_nu),e1out(3,3,0:n_nu)
      integer :: i,j,ii,jj
      real*8 :: pi8,e00m1,ee,sqpi8,cc,t1h(3,3),tl(3,3),tli(3,3),te(3,3)
      real*8, allocatable :: pvsq(:,:),avq(:,:),x21(:,:)
      allocate(pvsq(nne,3))
      allocate(avq(nne,3))
      do i=1,nne
        pvsq(i,:)=p_wing(:,i,ind_nu)*vsq(i)
      enddo
      avq=0.d0
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          do ii=1,3
            cc=pvsq(j1,ii)
            do i=1,nrow
              i1=refrow(i)
              avq(i1,ii)=avq(i1,ii)+w_s(i,j)*cc
            enddo
          enddo
        enddo
        call dgop(avq,3*nne,'  +',comm_b)
#endif
      else
        do j=1,nne
          do ii=1,3
            cc=pvsq(j,ii)
            do i=1,nne
              avq(i,ii)=avq(i,ii)+w_s(i,j)*cc
            enddo
          enddo
        enddo
      endif
c ------------------------------------------------------------------
      pi8=8.d0*pi
      sqpi8=sqrt(pi8)
      p00=0.d0
      do i=1,3
        p00=p00+p_head(i,i,ind_nu)
      enddo
      p00=p00/3.d0
c ------- Tensor 1-8*pi*H ------------------------------------------
      t1h=-pi8*p_head(:,:,ind_nu)
      do i=1,3
        t1h(i,i)=1.d0+t1h(i,i)
      enddo
c ------- Tensor L -------------------------------------------------
      tl=0.d0
      do jj=1,3
        do ii=1,3
          do i=1,nne
            tl(ii,jj)=tl(ii,jj)+pvsq(i,ii)*avq(i,jj)
          enddo
        enddo
      enddo
      tl=t1h-pi8*tl
c ------------------- 1/(e*L*e) ------------------------------------
      tli=tl
      call invers_r(3,tli,3)
      e00m1=0.d0
      do i=1,3
        e00m1=e00m1+tli(i,i)
      enddo
      e00m1=e00m1/3.d0
      e1out(:,:,i_nu)=tli
      te=tli
      do i=1,3
        te(i,i)=te(i,i)-1.d0
      enddo
      ws_head_tens_r(:,:,ind_nu)=pi8*te
c ------------------------------------------------------------------
      eout(:,:,i_nu)=tl
      ee=0.d0
      do i=1,3
        ee=ee+tl(i,i)
      enddo
      e_output(i_nu,1)=ee/3.d0
c ------------------------------------------------------------------
      call vec_tens_prod_r(tli,te)
      allocate(x21(nne,3))
      x21=0.d0
      do jj=1,3
        do ii=1,3
          x21(:,jj)=x21(:,jj)+avq(:,ii)*te(ii,jj)
        enddo
      enddo
      x21=sqpi8*x21
      do i=1,nne
        do ii=1,3
          ws_wing_r(ii,i,ind_nu)=x21(i,ii)*vsq(i)
        enddo
      enddo
      ws_wing_r(:,:,ind_nu)=sqpi8*ws_wing_r(:,:,ind_nu)
c ------------------------------------------------------------------
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          do i=1,nrow
            i1=refrow(i)
            cc=0.d0
            do ii=1,3
              cc=cc+avq(i1,ii)*x21(j1,ii)
            enddo
            cc=sqpi8*cc/3.d0
            w_s(i,j)=w_s(i,j)+cc
          enddo
        enddo
        do j=1,ncol
          j1=refcol(j)
          do i=1,nrow
            i1=refrow(i)
            if(i1==j1) w_s(i,j)=w_s(i,j)-1.d0
          enddo
        enddo
        do j=1,ncol
          j1=refcol(j)
          do i=1,nrow
            i1=refrow(i)
            w_s(i,j)=vsq(i1)*w_s(i,j)*vsq(j1)
          enddo
        enddo
#endif
      else
        do j=1,nne
          do i=1,nne
            cc=0.d0
            do ii=1,3
              cc=cc+avq(i,ii)*x21(j,ii)
            enddo
            cc=sqpi8*cc/3.d0
            w_s(i,j)=w_s(i,j)+cc
          enddo
        enddo
        do i=1,nne
          w_s(i,i)=w_s(i,i)-1.d0
        enddo
        do j=1,nne
          do i=1,nne
            w_s(i,j)=vsq(i)*w_s(i,j)*vsq(j)
          enddo
        enddo
      endif
      ws_head_nu(i_nu)=pi8*(e00m1-1.d0)
      deallocate(pvsq,avq,x21)
      end
