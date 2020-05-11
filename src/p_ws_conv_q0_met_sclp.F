      subroutine p_ws_conv_q0_met_sclp(n,nrow,ncol,pwing,wwing,p00,w00,
     &                                 p11,w11,pws,q2a)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol
      real*8, intent(in) :: p00,w00,q2a
      complex*16, intent(in) :: pwing(n),wwing(n),p11(nrow,ncol),
     &                          w11(nrow,ncol)
      real*8, intent(out) :: pws
      integer :: j,i
      real*8 :: d
      pws=q2a*p00*w00
      do i=1,n
        pws=pws+pwing(i)*conjg(wwing(i))
      enddo
      do i=1,n
        pws=pws+conjg(pwing(i))*wwing(i)
      enddo
      d=0.d0
      do i=1,nrow
        do j=1,ncol
          d=d+p11(i,j)*conjg(w11(i,j))
        enddo
      enddo
      call dgop(d,1,'  +',comm_b)
      pws=pws+d
      end


      subroutine p_ws_conv_q0_met_sclp_r(n,nrow,ncol,pwing,wwing,p00,
     &                                   w00,p11,w11,pws,q2a)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,nrow,ncol
      real*8, intent(in) :: p00,w00,q2a,pwing(n),wwing(n),
     &                      p11(nrow,ncol),w11(nrow,ncol)
      real*8, intent(out) :: pws
      integer :: j,i
      real*8 :: d
      pws=q2a*p00*w00
      do i=1,n
        pws=pws+pwing(i)*wwing(i)
      enddo
      do i=1,n
        pws=pws+pwing(i)*wwing(i)
      enddo
      d=0.d0
      do i=1,nrow
        do j=1,ncol
          d=d+p11(i,j)*w11(i,j)
        enddo
      enddo
      call dgop(d,1,'  +',comm_b)
      pws=pws+d
      end
