      subroutine wsq0_met(ncol,nrow,refcol,refrow,nblock,desc,nne,scal,
     &                    p01,vsq,p0,w_s,n1,n2)
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
     &                       nblock,desc(9),n1,n2
      real*8, intent(in) :: p0,vsq(nne)
      complex*16, intent(in) :: p01(nne)
      complex*16, intent(inout) :: w_s(n1,n2)
      integer :: i,j
      real*8 :: dd
      complex*16 :: cc
      complex*16, allocatable :: pvsq(:),apv(:)
      allocate(pvsq(nne))
      allocate(apv(nne))
      do i=1,nne
        pvsq(i)=p01(i)*vsq(i)
      enddo
      apv=(0.d0,0.d0)
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          cc=conjg(pvsq(j1))
          do i=1,nrow
            i1=refrow(i)
            apv(i1)=apv(i1)+w_s(i,j)*cc
          enddo
        enddo
        call dgop(apv,2*nne,'  +',comm_b)
#endif
      else
        do j=1,nne
          cc=conjg(pvsq(j))
          do i=1,nne
            apv(i)=apv(i)+w_s(i,j)*cc
          enddo
        enddo
      endif
      dd=p0
      do i=1,nne
        dd=dd+pvsq(i)*apv(i)
      enddo
      do i=1,nne
        ws01_met(i)=-conjg(apv(i))*vsq(i)/dd
      enddo
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          cc=conjg(apv(j1))/dd
          do i=1,nrow
            i1=refrow(i)
            w_s(i,j)=w_s(i,j)-apv(i1)*cc
          enddo
        enddo
#endif
      else
        do j=1,nne
          cc=conjg(apv(j))/dd
          do i=1,nne
            w_s(i,j)=w_s(i,j)-apv(i)*cc
          enddo
        enddo
      endif
      deallocate(pvsq,apv)
      end


      subroutine wsq0_met_r(ncol,nrow,refcol,refrow,nblock,desc,nne,
     &                      scal,p01,vsq,p0,w_s,n1,n2)
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
     &                       nblock,desc(9),n1,n2
      real*8, intent(in) :: p0,vsq(nne),p01(nne)
      real*8, intent(inout) :: w_s(n1,n2)
      integer :: i,j
      real*8 :: dd,cc
      real*8, allocatable :: pvsq(:),apv(:)
      allocate(pvsq(nne))
      allocate(apv(nne))
      do i=1,nne
        pvsq(i)=p01(i)*vsq(i)
      enddo
      apv=0.d0
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          cc=pvsq(j1)
          do i=1,nrow
            i1=refrow(i)
            apv(i1)=apv(i1)+w_s(i,j)*cc
          enddo
        enddo
        call dgop(apv,nne,'  +',comm_b)
#endif
      else
        do j=1,nne
          cc=pvsq(j)
          do i=1,nne
            apv(i)=apv(i)+w_s(i,j)*cc
          enddo
        enddo
      endif
      dd=p0
      do i=1,nne
        dd=dd+pvsq(i)*apv(i)
      enddo
      do i=1,nne
        ws01_met_r(i)=-apv(i)*vsq(i)/dd
      enddo
      if(scal) then
#ifdef SCLP
        do j=1,ncol
          j1=refcol(j)
          cc=apv(j1)/dd
          do i=1,nrow
            i1=refrow(i)
            w_s(i,j)=w_s(i,j)-apv(i1)*cc
          enddo
        enddo
#endif
      else
        do j=1,nne
          cc=apv(j)/dd
          do i=1,nne
            w_s(i,j)=w_s(i,j)-apv(i)*cc
          enddo
        enddo
      endif
      deallocate(pvsq,apv)
      end




      subroutine wsq0_met_rO(nne,p01,vsq,p0,w_s)
      use atom_mod
      use etot_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nne
      real*8, intent(in) :: p0,vsq(nne),p01(nne)
      real*8, intent(inout) :: w_s(nne,nne)
      integer :: i,j
      real*8 :: dd,cc
      real*8, allocatable :: pvsq(:),apv(:)
      allocate(pvsq(nne))
      allocate(apv(nne))
      do i=1,nne
        pvsq(i)=p01(i)*vsq(i)
      enddo
      apv=0.d0
      do j=1,nne
        cc=pvsq(j)
        do i=1,nne
          apv(i)=apv(i)+w_s(i,j)*cc
        enddo
      enddo
      dd=p0
      do i=1,nne
        dd=dd+pvsq(i)*apv(i)
      enddo
      do i=1,nne
        ws01_met_r(i)=-apv(i)*vsq(i)/dd
      enddo
      do j=1,nne
        cc=apv(j)/dd
        do i=1,nne
          w_s(i,j)=w_s(i,j)-apv(i)*cc
        enddo
      enddo
      deallocate(pvsq,apv)
      end
