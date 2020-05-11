      subroutine interp_g_wmt(gexa,gspl_geom,gspl_asy,gasy,f_out,w1,
     &                        n,n1,indi,indj,key)
c ------- We assume that G(-w) = G*(w), i.e. G(t) - real -----
c     key = 1 - for  i a/w + b/w^2
c     key = 3 - for  i a/w^3 + b/w^4
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,key,n1,indi,indj
      real*8, intent(in) :: gexa(nfun_red,nfun_red,2,0:n_omega_exa),
     &                      gspl_geom(nfun_red,nfun_red,2,
     &                                0:n_omega_geom+2,4),w1,
     &                      gspl_asy(nfun_red,nfun_red,2,n_omega_asy+1,
     &                               4),
     &                      gasy(nfun_red,nfun_red,2)
      complex*16, intent(out) :: f_out(n,n1)
      integer :: i,j,i0,iw,ii
      real*8 :: w,a,b,pib,x,x2,x3,xm
      real*8, allocatable :: ab(:,:,:)
      pib=pi/betta_t
      w=abs(w1)
      iw=(idnint(w/pib)-1)/2
      allocate(ab(n,n1,2))
      if(iw<=n_omega_exa) then
        ab=gexa(indi:indi+n-1,indj:indj+n1-1,:,iw)
      else if(iw<=mats_num_omega(n_omega_exa+n_omega_geom+1)) then
        ii=n_omega_geom+1
        do i0=n_omega_exa+1,n_omega_exa+n_omega_geom
          if(mats_num_omega(i0)>=iw) then
            ii=i0-n_omega_exa
            goto 1
          endif
        enddo
1       x=w-w_omega(n_omega_exa+ii-1)
        x2=x*x
        x3=x2*x
        do j=indj,indj+n1-1
          do i=indi,indi+n-1
            a=gspl_geom(i,j,1,ii,1)+gspl_geom(i,j,1,ii,2)*x
     &       +gspl_geom(i,j,1,ii,3)*x2+gspl_geom(i,j,1,ii,4)*x3
            b=gspl_geom(i,j,2,ii,1)+gspl_geom(i,j,2,ii,2)*x
     &       +gspl_geom(i,j,2,ii,3)*x2+gspl_geom(i,j,2,ii,4)*x3
            ab(i-indi+1,j-indj+1,1)=a
            ab(i-indi+1,j-indj+1,2)=b
          enddo
        enddo
      else if(iw<=mats_num_omega(n_omega-1)) then
        ii=n_omega_asy
        do i0=n_omega_exa+n_omega_geom+2,n_omega-2
          if(mats_num_omega(i0)>=iw) then
            ii=i0-n_omega_exa-n_omega_geom
            goto 2
          endif
        enddo
2       xm=(1.d0/w**2-bb_omega)/aa_omega
        x=xm-xm_omega(ii-1)
        x2=x*x
        x3=x2*x
        do j=indj,indj+n1-1
          do i=indi,indi+n-1
            a=gspl_asy(i,j,1,ii,1)+gspl_asy(i,j,1,ii,2)*x
     &       +gspl_asy(i,j,1,ii,3)*x2+gspl_asy(i,j,1,ii,4)*x3
            b=gspl_asy(i,j,2,ii,1)+gspl_asy(i,j,2,ii,2)*x
     &       +gspl_asy(i,j,2,ii,3)*x2+gspl_asy(i,j,2,ii,4)*x3
            ab(i-indi+1,j-indj+1,1)=a
            ab(i-indi+1,j-indj+1,2)=b
          enddo
        enddo
      else
        if(key==1) then
          do j=indj,indj+n1-1
            do i=indi,indi+n-1
              a=gasy(i,j,1)*w_omega(n_omega)**2/w**2
              b=gasy(i,j,2)*w_omega(n_omega)/w
              ab(i-indi+1,j-indj+1,1)=a
              ab(i-indi+1,j-indj+1,2)=b
            enddo
          enddo
        else if(key==3) then
          do j=indj,indj+n1-1
            do i=indi,indi+n-1
              a=gasy(i,j,1)*w_omega(n_omega)**4/w**4
              b=gasy(i,j,2)*w_omega(n_omega)**3/w**3
              ab(i-indi+1,j-indj+1,1)=a
              ab(i-indi+1,j-indj+1,2)=b
            enddo
          enddo
        endif
      endif
      f_out=dcmplx(ab(:,:,1),ab(:,:,2))
      deallocate(ab)
      end
