      subroutine interp_g_omega_sq(f_in,f_out,x_out,n_out,n0,n1,n)
c ------- We assume that G_ij(-w) = G_ji*(w), i.e. G(t) - hermitian -----
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n0,n1,n
      real*8, intent(in) :: x_out(0:n_out)
	complex*16, intent(in) :: f_in(n0,n0,0:n_omega)
	complex*16, intent(out) :: f_out(n1,n1,0:n_out)
	integer :: i,io,in
	real*8 :: v,v2,v4,p2b,pb,x,x_cheb_from_omega
	real*8, allocatable :: tmp(:,:,:,:),chb(:,:,:,:),cheb(:),ab(:,:,:)
	f_out=(0.d0,0.d0)
	pb=pi/betta_t
	p2b=pb+pb
	do io=0,n_out
	  v=x_out(io)
	  if(x_out(io)<0) v=-v
	  if(v<=w_omega(n_omega_exa)) then
	    in=idnint((v-pb)/p2b)
	    v2=pb+in*p2b
	    if(abs(v-v2)<1.d-13) then
	      f_out(1:n,1:n,io)=f_in(1:n,1:n,in)
	    else
	      call interp_array_opt(f_in,f_out(1,1,io),w_omega,v,
     &	                          n_omega_exa,0,2*n0,n0,2*n1,n1,2*n,n,5)
          endif
	  else
	    allocate(tmp(n,n,2,n_omega_asy))
	    allocate(ab(n,n,2))
	    do i=1,n_omega_asy
	      call ferm_pack_omega(f_in(1,1,n_omega_exa+i),tmp(1,1,1,i),n,
     &	                         n0,n)
          enddo
	    if(v<=w_omega(n_omega)) then
	      allocate(chb(n,n,2,0:n_omega_asy-2))
	      call dgemm('n','t',2*n*n,n_omega_asy-1,n_omega_asy-1,
     &	               1.d0,tmp,2*n*n,cheb_omega(0,1),n_omega_asy-1,
     &                 0.d0,chb,2*n*n)
	      allocate(cheb(0:n_omega_asy-2))
            x=x_cheb_from_omega(v)
	      call dcheb_recurr(x,cheb,n_omega_asy-2)
	      ab=0.d0
	      do i=0,n_omega_asy-2
	        ab=ab+chb(:,:,:,i)*cheb(i)
	      enddo
	      deallocate(cheb,chb)
	    else
            v2=v*v
            v4=v2*v2
	      ab(:,:,1)=(asy_omega_interp(1,1,1)*tmp(:,:,1,n_omega_asy-1)
     &	             +asy_omega_interp(2,1,1)*tmp(:,:,1,n_omega_asy))/v2
     &               +(asy_omega_interp(1,2,1)*tmp(:,:,1,n_omega_asy-1)
     &	             +asy_omega_interp(2,2,1)*tmp(:,:,1,n_omega_asy))/v4
            v2=v
            v4=v*v2
	      ab(:,:,2)=(asy_omega_interp(1,1,2)*tmp(:,:,2,n_omega_asy-1)
     &	             +asy_omega_interp(2,1,2)*tmp(:,:,2,n_omega_asy))/v2
     &               +(asy_omega_interp(1,2,2)*tmp(:,:,2,n_omega_asy-1)
     &	             +asy_omega_interp(2,2,2)*tmp(:,:,2,n_omega_asy))/v4
	    endif
	    call ferm_unpack_omega(f_out(1,1,io),ab,n,n,n1)
	    deallocate(tmp,ab)
	  endif
	  if(x_out(io)<0) f_out(:,:,io)=conjg(transpose(f_out(:,:,io)))
	enddo   !! over io
      end




      subroutine interp_g_sq_spl(gexa,gspl_geom,gspl_asy,gasy,f_out,w1,
     &                           n,key)
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
      real*8, intent(in) :: gexa(n,n,2,0:n_omega_exa),
     &                      gspl_geom(n,n,2,0:n_omega_geom+2,4),w1,
     &                      gspl_asy(n,n,2,n_omega_asy+1,4),gasy(n,n,2)
      complex*16, intent(out) :: f_out(n,n)
      integer :: i,j,i0,iw,ii
      real*8 :: w,a,b,pib,x,x2,x3,xm
      real*8, allocatable :: ab(:,:,:)
      pib=pi/betta_t
      w=abs(w1)
      iw=(idnint(w/pib)-1)/2
      allocate(ab(n,n,2))
      if(iw<=n_omega_exa) then
        ab=gexa(:,:,:,iw)
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
        do j=1,n
          do i=1,n
            a=gspl_geom(i,j,1,ii,1)+gspl_geom(i,j,1,ii,2)*x
     &       +gspl_geom(i,j,1,ii,3)*x2+gspl_geom(i,j,1,ii,4)*x3
            b=gspl_geom(i,j,2,ii,1)+gspl_geom(i,j,2,ii,2)*x
     &       +gspl_geom(i,j,2,ii,3)*x2+gspl_geom(i,j,2,ii,4)*x3
            ab(i,j,1)=a
            ab(i,j,2)=b
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
        do j=1,n
          do i=1,n
            a=gspl_asy(i,j,1,ii,1)+gspl_asy(i,j,1,ii,2)*x
     &       +gspl_asy(i,j,1,ii,3)*x2+gspl_asy(i,j,1,ii,4)*x3
            b=gspl_asy(i,j,2,ii,1)+gspl_asy(i,j,2,ii,2)*x
     &       +gspl_asy(i,j,2,ii,3)*x2+gspl_asy(i,j,2,ii,4)*x3
            ab(i,j,1)=a
            ab(i,j,2)=b
          enddo
        enddo
      else
        if(key==1) then
          do j=1,n
            do i=1,n
              a=gasy(i,j,1)*w_omega(n_omega)**2/w**2
              b=gasy(i,j,2)*w_omega(n_omega)/w
              ab(i,j,1)=a
              ab(i,j,2)=b
            enddo
          enddo
        else if(key==3) then
          do j=1,n
            do i=1,n
              a=gasy(i,j,1)*w_omega(n_omega)**4/w**4
              b=gasy(i,j,2)*w_omega(n_omega)**3/w**3
              ab(i,j,1)=a
              ab(i,j,2)=b
            enddo
          enddo
        endif
      endif
	call ferm_unpack_omega(f_out,ab,n,n,n)
      deallocate(ab)
      if(w1<0) f_out(:,:)=conjg(transpose(f_out(:,:)))
      end