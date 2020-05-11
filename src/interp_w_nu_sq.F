      subroutine interp_w_nu_sq_spl(wexa,wspl_geom,wspl_asy,wasy,f_out,
     &                              w1,n)
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n
      real*8, intent(in) :: wexa(n,n,0:n_nu_exa),
     &                      wspl_geom(n,n,0:n_nu_geom+2,4),w1,
     &                      wspl_asy(n,n,n_nu_asy+1,4),wasy(n,n)
      complex*16, intent(out) :: f_out(n,n)
      real*8, allocatable :: tmp(:,:)
      integer :: i,i0,iw
      real*8 :: w,pib,x,x2,x3,xm
      pib=pi/betta_t
      w=abs(w1)
      iw=idnint(w/pib/2)
      allocate(tmp(n,n))
      if(iw<=n_nu_exa) then
        tmp=wexa(:,:,iw)
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
        tmp=wspl_geom(:,:,i,1)+wspl_geom(:,:,i,2)*x
     &     +wspl_geom(:,:,i,3)*x2+wspl_geom(:,:,i,4)*x3
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
        tmp=wspl_asy(:,:,i,1)+wspl_asy(:,:,i,2)*x
     &     +wspl_asy(:,:,i,3)*x2+wspl_asy(:,:,i,4)*x3
      else
        tmp=wasy(:,:)*w_nu(n_nu)**2/w**2
      endif
	  call unpack_hermit(f_out,tmp,n,n,n,(0.d0,0.d0),(1.d0,0.d0))
      deallocate(tmp)
      end
      
      
      
      subroutine interp_w_nu_re_spl(wexa,wspl_geom,wspl_asy,wasy,f_out,
     &                              w1,n,n0)
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0
      real*8, intent(in) :: wexa(n0,n0,0:n_nu_exa),
     &                      wspl_geom(n0,n0,0:n_nu_geom+2,4),w1,
     &                      wspl_asy(n0,n0,n_nu_asy+1,4),wasy(n0,n0)
      complex*16, intent(out) :: f_out(n0,n0)
      integer :: i,i0,iw
      real*8 :: w,pib,x,x2,x3,xm
      real*8, allocatable :: tmp(:,:)
      pib=pi/betta_t
      w=abs(w1)
      iw=idnint(w/pib/2)
      allocate(tmp(n,n))
      if(iw<=n_nu_exa) then
        tmp=wexa(1:n,1:n,iw)
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
        tmp=wspl_geom(1:n,1:n,i,1)+wspl_geom(1:n,1:n,i,2)*x
     &     +wspl_geom(1:n,1:n,i,3)*x2+wspl_geom(1:n,1:n,i,4)*x3
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
        tmp=wspl_asy(1:n,1:n,i,1)+wspl_asy(1:n,1:n,i,2)*x
     &     +wspl_asy(1:n,1:n,i,3)*x2+wspl_asy(1:n,1:n,i,4)*x3
      else
        tmp=wasy(1:n,1:n)*w_nu(n_nu)**2/w**2
      endif
	  f_out(1:n,1:n)=tmp
      deallocate(tmp)
      end
