      subroutine decode_generators(t,nt,g,g_j,a,ifiaa,lmax,jmax,ng,h,
     &                             h_j,su1,su1_j,fac)
c     Parses the string t containing the group element describers
c     and makes a matrix g and vector a for each element.
c     space group ops are defined as:  (g,a)(p)=g*p+a
c
c     Point group element describers:
c     Rnv = n-fold rot around vec v, n pos integer   ) also products
c     Mv  = mirror taking vec v into -v              ) of form 'a*b' of
c     I   = inversion                                ) any two
c     -   = J => -J                                  ) of these
c
c     Glide planes, screw axes etc:
c     append  ':Tv'  to any point operation above where v is the
c     translation vector for the operation.
c
c     Vectors v are all given as string-vector (x,y,z)
c           (see formats below)
c     or abbreviate   'D' for (1,1,1), 'X' for (1,0,0),  Y,Z  similar.
c     Possible formats for numbers x,y,z:
c                  a
c                  a/b
c                  sqrt(a)
c                  sqrt(a/b)
c              with a,b integer or real.
      use manager_mod
      use parallel_mod
      use solid_mod
      character*80, intent(in) :: t
      integer, intent(in) :: nt,lmax,jmax
      real*8, intent(in) :: fac(0:*)
      integer, intent(out) :: ifiaa(*),ng
      real*8, intent(out) :: g(maxwig,*),a(3,*),h(maxwig),
     &                       su1(-lmax:lmax,-lmax:lmax)
      complex*16, intent(out) :: g_j(maxwig_j,*),h_j(maxwig_j),
     &                           su1_j(jmax+1,jmax+1)
      integer :: i,l,iwig,ist,ii,jj
      ng=0
      i=1
  90  call parsbl(t,nt,i)
      if(i.gt.nt) return
      ng=ng+1
      ifiaa(ng)=1
      if(t(i:i).eq.'-') then
        ifiaa(ng)=-1
        i=i+1
      endif
      call decode_point_op(t,nt,i,g(1,ng),g_j(1,ng),lmax,jmax,su1,su1_j,
     &                     fac)
      if(t(i:i)=='*') then
        i=i+1
        call decode_point_op(t,nt,i,h,h_j,lmax,jmax,su1,su1_j,fac)
        do l=0,lmax
          iwig=l*(2*l-1)*(2*l+1)/3+1
          call dcopy((2*l+1)**2,g(iwig,ng),1,su1(-lmax,-lmax),1)
          call prd(su1(-lmax,-lmax),h(iwig),g(iwig,ng),l)
        enddo
        if(irel.eq.2) then
          do l=0,lmax
            ist=-1
            if(l.eq.0) ist=1
            do ii=ist,1,2
              jj=l+l+ii
              iwig=4*l*(2*l*l+1)/3+2*(ii+1)*l*l+1
              call zcopy((jj+1)**2,g_j(iwig,ng),1,su1_j,1)
              call prd_j(su1_j,h_j(iwig),g_j(iwig,ng),jj)
            enddo
          enddo
        endif
      endif
      call dcopy(3,0.d0,0,a(1,ng),1)
      if(t(i:i).eq.':') then
        i=i+1
        if(t(i:i).ne.'T') then
          if(maswrk) write(iun,*)'*** parsgn: "T" expected after ":"'
          call ending
        endif
        i=i+1
        call parsvc(t,nt,i,a(1,ng),rbas)
      endif
      if(t(i:i).ne.'_') then
        if(maswrk) write(iun,*)'*** no blank found after op'
        call ending
      endif
      goto 90
      end
