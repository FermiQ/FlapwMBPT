      subroutine var_int_b(ww,p,q,p0,q0,f,g,fg,r,bl,cl,work,nr,h)
      use manager_mod
      use units_mod
      implicit none
      integer :: nr,ir
      real*8 ::  ww(nr),p(nr),q(nr),p0(nr),q0(nr),f(nr),g(nr),fg(nr),
     &           bl(nr),cl(nr),r(nr),work(nr),c2,h
      c2=clight**2
      do ir=1,nr
        work(ir)=p(ir)*p0(ir)
      enddo
      if(irel>=1) then
        do ir=1,nr
          work(ir)=work(ir)+q(ir)*q0(ir)/c2
        enddo
      endif
      do ir=1,nr
        work(ir)=ww(ir)*work(ir)*h
      enddo
      do ir=1,nr
        f(ir)=work(ir)*bl(ir)
      enddo
      call intalt(nr,f,g)
      fg(1)=0.d0
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
