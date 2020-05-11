      subroutine interp_w_nu_spl(wexa,wspl_geom,wspl_asy,wasy,f_out,w1,
     &                           n)
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: wexa(0:n_nu_exa,n),
     &                      wspl_geom(0:n_nu_geom+2,4,n),w1,
     &                      wspl_asy(n_nu_asy+1,4,n),wasy(n)
      real*8, intent(out) :: f_out(n)
      integer :: i,j,i0,iw
      real*8 :: w,pib,x,x2,x3,xm
      pib=pi/betta_t
      w=abs(w1)
      iw=idnint(w/pib/2)
      if(iw<=n_nu_exa) then
        f_out=wexa(iw,:)
      else if(iw<=mats_num_nu(n_nu_exa+n_nu_geom+1)) then
        i=n_nu_geom+1
        do i0=n_nu_exa+1,n_nu_exa+n_nu_geom
          if(mats_num_nu(i0)>=iw) then
            i=i0-n_nu_exa
            goto 1
          endif
        enddo
1       x=w-w_nu(n_nu_exa+i-1)
        x2=x*x
        x3=x2*x
        do j=1,n
          f_out(j)=wspl_geom(i,1,j)+wspl_geom(i,2,j)*x
     &            +wspl_geom(i,3,j)*x2+wspl_geom(i,4,j)*x3
        enddo
      else if(iw<=mats_num_nu(n_nu-1)) then
        i=n_nu_asy
        do i0=n_nu_exa+n_nu_geom+2,n_nu-2
          if(mats_num_nu(i0)>=iw) then
            i=i0-n_nu_exa-n_nu_geom
            goto 2
          endif
        enddo
2       xm=(1.d0/w**2-bb_nu)/aa_nu
        x=xm-xm_nu(i-1)
        x2=x*x
        x3=x2*x
        do j=1,n
          f_out(j)=wspl_asy(i,1,j)+wspl_asy(i,2,j)*x
     &            +wspl_asy(i,3,j)*x2+wspl_asy(i,4,j)*x3
        enddo
      else
        do j=1,n
          f_out(j)=wasy(j)*w_nu(n_nu)**2/w**2
        enddo
      endif
      end