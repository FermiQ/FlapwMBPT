      subroutine var_int1_b(ww,p,f,g,fg,r,bl,cl,work,nr,h)
      implicit none
      integer, intent(in) :: nr
      real*8, intent(in) :: h,ww(nr),p(nr),r(nr),bl(nr),cl(nr)
      real*8, intent(inout) :: f(nr),g(nr),work(nr)
      real*8, intent(out) :: fg(nr)
      integer :: ir
      do ir=1,nr
        work(ir)=ww(ir)*p(ir)*h*r(ir)**2
      enddo
      do ir=1,nr
        f(ir)=work(ir)*bl(ir)
      enddo
      call intalt(nr,f,g)
      do ir=2,nr
        fg(ir)=g(ir)*cl(ir)
      enddo
      f(1)=0.d0
      do ir=2,nr
        f(ir)=work(ir)*cl(ir)
      enddo
      call intalt(nr,f,g)
      do ir=1,nr
        fg(ir)=fg(ir)+(g(nr)-g(ir))*bl(ir)
      enddo
      end
