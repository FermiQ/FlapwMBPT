      subroutine decode_number(t,i1,i2,r)
      character*80, intent(in) :: t
      integer, intent(in) :: i1,i2
      real*8, intent(out) :: r
      integer :: j1,j2,n,isig,j,m
      real*8 :: r1
      j=i1
      isig=1
      if(t(i1:i1)=='-') then
        j=i1+1
        isig=-1
      endif
      if(t(j:j)=='s') then
        j1=j+5
        j2=i2-1
        m=index(t(j1:j2),'/')
        n=j1+m-1
        if(m==0) then
          read(t(j1:j2),*)r
          r=sqrt(r)
        else
          read(t(j1:n-1),*)r
          read(t(n+1:j2),*)r1
          r=sqrt(r/r1)
        endif
      else
        j1=j
        j2=i2
        m=index(t(j1:j2),'/')
        n=j1+m-1
        if(m==0) then
          read(t(j1:j2),*)r
        else
          read(t(j1:n-1),*)r
          read(t(n+1:j2),*)r1
          r=r/r1
        endif
      endif
      r=isig*r
      end
