      subroutine interp_p_nu_spl(pexa,pspl_geom,pspl_asy,pasy,f_out,w1,
     &                           n)
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: pexa(0:n_nu_exa,n),
     &                      pspl_geom(0:n_nu_geom+2,4,n),w1,
     &                      pspl_asy(n_nu_asy+1,4,n),pasy(n)
      real*8, intent(out) :: f_out(n)
      integer :: i,j,i0,iw
      real*8 :: w,pib,x,x2,x3,xm
      pib=pi/betta_t
      w=abs(w1)
      iw=idnint(w/pib/2)
      if(iw<=n_nu_exa) then
        f_out=pexa(iw,:)
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
          f_out(j)=pspl_geom(i,1,j)+pspl_geom(i,2,j)*x
     &            +pspl_geom(i,3,j)*x2+pspl_geom(i,4,j)*x3
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
          f_out(j)=pspl_asy(i,1,j)+pspl_asy(i,2,j)*x
     &            +pspl_asy(i,3,j)*x2+pspl_asy(i,4,j)*x3
        enddo
      else
        do j=1,n
          f_out(j)=pasy(j)*w_nu(n_nu)**2/w**2
        enddo
      endif
      end
      
      
      
      
      subroutine interp_p_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,f_out,
     &                              w1,n,n0)
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0
      real*8, intent(in) :: wexa(n,n,0:n_nu_exa),
     &                      wspl_geom(n,n,0:n_nu_geom+2,4),w1,
     &                      wspl_asy(n,n,n_nu_asy+1,4),wasy(n,n)
      real*8, intent(out) :: f_out(n0,n0)
      integer :: i,i0,iw
      real*8 :: w,pib,x,x2,x3,xm
      pib=pi/betta_t
      w=abs(w1)
      iw=idnint(w/pib/2)
      if(iw<=n_nu_exa) then
        f_out(1:n,1:n)=wexa(:,:,iw)
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
        f_out(1:n,1:n)=wspl_geom(:,:,i,1)+wspl_geom(:,:,i,2)*x
     &                +wspl_geom(:,:,i,3)*x2+wspl_geom(:,:,i,4)*x3
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
        f_out(1:n,1:n)=wspl_asy(:,:,i,1)+wspl_asy(:,:,i,2)*x
     &                +wspl_asy(:,:,i,3)*x2+wspl_asy(:,:,i,4)*x3
      else
        f_out(1:n,1:n)=wasy(:,:)*w_nu(n_nu)**2/w**2
      endif
      end