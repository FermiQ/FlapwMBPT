	  subroutine coarse_interp_ac
	  use manager_mod
	  use models_mod
	  use parallel_mod
	  use solid_mod
	  use vertex_mod
	  implicit none
	  integer :: k,i,j,n,n0,nnn,n5,i1,i2,i3
	  real*8 :: v(3),d0,dd
	  integer, allocatable :: ifl(:)
	  real*8, allocatable :: pnt0(:,:),dist(:),r8(:)
c	  We define barycentric coordinates for each NPNT ------------
      nnn=maxval(ndiv_c)*2
      n5=(2*nnn+1)**3
      num_baryc_ac_max=20
      n0=4
	  allocate(pnt_baryc_ac(3,num_baryc_ac_max,npnt))
	  allocate(coef_baryc_ac(num_baryc_ac_max,npnt))
	  allocate(num_baryc_ac(npnt))
	  allocate(ifl(n5))
      allocate(pnt0(3,n5))
      allocate(dist(n5))
      allocate(r8(num_baryc_ac_max))
	  do k=1,npnt
c --- First we find closiest to K point on coarse mesh ---------	
	    j=0
        do i3=-nnn,nnn
          do i2=-nnn,nnn
            do i1=-nnn,nnn
	          j=j+1
	          pnt0(:,j)=i1*qb0_c(:,1)+i2*qb0_c(:,2)+i3*qb0_c(:,3)
	          v=pnt0(:,j)-pnt(:,k)
	          dist(j)=dot_product(v,v)
	        enddo
	      enddo
	    enddo
	    call ordering_all(n5,dist,ifl)
c ----------------------------------------------------------------
        if(dist(ifl(1))<1.d-8) then
          num_baryc_ac(k)=1
          coef_baryc_ac(1,k)=1.d0
          pnt_baryc_ac(:,1,k)=pnt0(:,ifl(1))
        else
          d0=dist(ifl(n0))
          n=n0
          do i=n0+1,n5
            if(abs(dist(ifl(i))-d0)<1.d-8) n=n+1
          enddo
          num_baryc_ac(k)=n
          dd=0.d0
          do i=1,n
            r8(i)=dist(ifl(i))
            dd=dd+1.d0/r8(i)
          enddo
          d0=1.d0/dd
          do i=1,n
            coef_baryc_ac(i,k)=d0/r8(i)
            pnt_baryc_ac(:,i,k)=pnt0(:,ifl(i))
          enddo
        endif
	  enddo    !! over k
	  deallocate(ifl,pnt0,dist,r8)
      num_baryc_ac_max=maxval(num_baryc_ac)
      allocate(ind_kk_k_ac(num_baryc_ac_max,npnt))
      allocate(pnt0(3,npnt*num_baryc_ac_max))
      num_kk_k_ac=0
      do k=1,npnt
        do i=1,num_baryc_ac(k)
          v=pnt_baryc_ac(:,i,k)-pnt(:,k)
          i2=0
          do i1=1,num_kk_k_ac
            dd=maxval(abs(pnt0(:,i1)-v))
            if(dd<1.d-8) i2=i1
          enddo
          if(i2==0) then
            num_kk_k_ac=num_kk_k_ac+1
            ind_kk_k_ac(i,k)=num_kk_k_ac
            pnt0(:,num_kk_k_ac)=v
          else
            ind_kk_k_ac(i,k)=i2
          endif
        enddo
      enddo
      allocate(dif_kk_k_ac(3,num_kk_k_ac))
      dif_kk_k_ac=pnt0(:,1:num_kk_k_ac)
      deallocate(pnt0)
      end
