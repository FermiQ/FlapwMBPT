      subroutine p_ws_conv_q0_sclp(n,nrow,ncol,pwing,wwing,p00,w00,p11,
     &                             w11,pws)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol
      complex*16, intent(in) :: p00(3,3),w00(3,3)
      complex*16, intent(in) :: pwing(3,n),wwing(3,n),p11(nrow,ncol),
     &                          w11(nrow,ncol)
      real*8, intent(out) :: pws
      integer :: j,i
      real*8 :: tt
      complex*16 :: t(3,3)
      call tensor_prod_via_ylm(p00,w00,t,4)
      pws=(t(1,1)+t(2,2)+t(3,3))/3
      tt=0.d0
      do j=1,n
        do i=1,3
          tt=tt+pwing(i,j)*conjg(wwing(i,j))
        enddo
      enddo
      pws=pws+2.d0*tt
      tt=0.d0
      do i=1,nrow
        do j=1,ncol
          tt=tt+p11(i,j)*conjg(w11(i,j))
        enddo
      enddo
      call dgop(tt,1,'  +',comm_b)
      pws=pws+tt
      end


      subroutine p_ws_conv_q0_sclp_r(n,nrow,ncol,pwing,wwing,p00,w00,
     &                               p11,w11,pws)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol
      real*8, intent(in) :: p00(3,3),w00(3,3),pwing(3,n),wwing(3,n),
     &                      p11(nrow,ncol),w11(nrow,ncol)
      real*8, intent(out) :: pws
      integer :: j,i
      real*8 :: tt
      real*8 :: t(3,3)
      call tensor_prod_via_ylm_r(p00,w00,t,4)
      pws=(t(1,1)+t(2,2)+t(3,3))/3
      tt=0.d0
      do j=1,n
        do i=1,3
          tt=tt+pwing(i,j)*wwing(i,j)
        enddo
      enddo
      pws=pws+2.d0*tt
      tt=0.d0
      do i=1,nrow
        do j=1,ncol
          tt=tt+p11(i,j)*w11(i,j)
        enddo
      enddo
      call dgop(tt,1,'  +',comm_b)
      pws=pws+tt
      end
