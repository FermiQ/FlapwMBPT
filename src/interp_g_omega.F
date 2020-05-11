      subroutine interp_g_omega(ab_in,ab_cheb,f_out,x_out,n_out,m,n,key)
c ------- We assume that G_ij(-w) = G_ij*(w), i.e. G(t) - real -----
c     key = 0 - for Gx+Gc
c     key = 1 - only for Gc
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n_out,n,m,key
      real*8, intent(in) :: ab_in(m,n,2,0:n_omega),x_out(0:n_out),
     &                      ab_cheb(m,n,2,0:n_omega_asy-2)
      complex*16, intent(out) :: f_out(m,n,0:n_out)
      integer :: i,j,io
      real*8 :: w,w1
      real*8, allocatable :: aa(:,:,:)
      allocate(aa(m,n,2))
      do io=0,n_out
      w1=x_out(io)
      w=abs(w1)
      call ferm_interp_omega(m*n,ab_in,ab_cheb,w,aa,key)
      do j=1,n
        do i=1,m
          f_out(i,j,io)=dcmplx(aa(i,j,1),aa(i,j,2))
        enddo
      enddo
      if(w1<0) f_out(:,:,io)=conjg(f_out(:,:,io))
      enddo   !! over io
      deallocate(aa)
      end




      subroutine interp_g_spl(gexa,gspl_geom,gspl_asy,gasy,f_out,w1,n,
     &                        key)
c ------- We assume that G(-w) = G*(w), i.e. G(t) - real -----
c     key = 1 - for  i a/w + b/w^2
c     key = 3 - for  i a/w^3 + b/w^4
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,key
      real*8, intent(in) :: gexa(0:n_omega_exa,2,n),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n),w1,
     &                      gspl_asy(n_omega_asy+1,4,2,n),gasy(2,n)
      complex*16, intent(out) :: f_out(n)
      integer :: i,j,i0,iw
      real*8 :: w,a,b,pib,x,x2,x3,xm
      pib=pi/betta_t
      w=abs(w1)
      iw=(idnint(w/pib)-1)/2
      if(iw<=n_omega_exa) then
        f_out=dcmplx(gexa(iw,1,:),gexa(iw,2,:))
      else if(iw<=mats_num_omega(n_omega_exa+n_omega_geom+1)) then
        i=n_omega_geom+1
        do i0=n_omega_exa+1,n_omega_exa+n_omega_geom
          if(mats_num_omega(i0)>=iw) then
            i=i0-n_omega_exa
            goto 1
          endif
        enddo
1       x=w-w_omega(n_omega_exa+i-1)
        x2=x*x
        x3=x2*x
        do j=1,n
          a=gspl_geom(i,1,1,j)+gspl_geom(i,2,1,j)*x
     &     +gspl_geom(i,3,1,j)*x2+gspl_geom(i,4,1,j)*x3
          b=gspl_geom(i,1,2,j)+gspl_geom(i,2,2,j)*x
     &     +gspl_geom(i,3,2,j)*x2+gspl_geom(i,4,2,j)*x3
          f_out(j)=dcmplx(a,b)
        enddo
      else if(iw<=mats_num_omega(n_omega-1)) then
        i=n_omega_asy
        do i0=n_omega_exa+n_omega_geom+2,n_omega-2
          if(mats_num_omega(i0)>=iw) then
            i=i0-n_omega_exa-n_omega_geom
            goto 2
          endif
        enddo
2       xm=(1.d0/w**2-bb_omega)/aa_omega
        x=xm-xm_omega(i-1)
        x2=x*x
        x3=x2*x
        do j=1,n
          a=gspl_asy(i,1,1,j)+gspl_asy(i,2,1,j)*x
     &     +gspl_asy(i,3,1,j)*x2+gspl_asy(i,4,1,j)*x3
          b=gspl_asy(i,1,2,j)+gspl_asy(i,2,2,j)*x
     &     +gspl_asy(i,3,2,j)*x2+gspl_asy(i,4,2,j)*x3
          f_out(j)=dcmplx(a,b)
        enddo
      else
        if(key==1) then
          do j=1,n
            a=gasy(1,j)*w_omega(n_omega)**2/w**2
            b=gasy(2,j)*w_omega(n_omega)/w
            f_out(j)=dcmplx(a,b)
          enddo
        else if(key==3) then
          do j=1,n
            a=gasy(1,j)*w_omega(n_omega)**4/w**4
            b=gasy(2,j)*w_omega(n_omega)**3/w**3
            f_out(j)=dcmplx(a,b)
          enddo
        endif
      endif
      if(w1<0) f_out(:)=conjg(f_out(:))
      end