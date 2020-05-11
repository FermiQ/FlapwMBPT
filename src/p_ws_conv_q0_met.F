      subroutine p_ws_conv_q0_met(n,n0,n1,pwing,wwing,p00,w00,p11,w11,
     &                            pws,q2a)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      real*8, intent(in) :: p00,w00,q2a
      complex*16, intent(in) :: pwing(n),wwing(n),p11(n0,n0),w11(n1,n1)
      real*8, intent(out) :: pws
      integer :: j,i
      pws=q2a*p00*w00
      do i=1,n
        pws=pws+pwing(i)*conjg(wwing(i))
      enddo
      do i=1,n
        pws=pws+conjg(pwing(i))*wwing(i)
        do j=1,n
          pws=pws+p11(i,j)*w11(j,i)
        enddo
      enddo
      end


      subroutine p_ws_conv_q0_met_r(n,n0,n1,pwing,wwing,p00,w00,p11,w11,
     &                              pws,q2a)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      real*8, intent(in) :: p00,w00,q2a,pwing(n),wwing(n),p11(n0,n0),
     &                      w11(n1,n1)
      real*8, intent(out) :: pws
      integer :: j,i
      pws=q2a*p00*w00
      do i=1,n
        pws=pws+pwing(i)*wwing(i)
      enddo
      do i=1,n
        pws=pws+pwing(i)*wwing(i)
        do j=1,n
          pws=pws+p11(i,j)*w11(j,i)
        enddo
      enddo
      end
