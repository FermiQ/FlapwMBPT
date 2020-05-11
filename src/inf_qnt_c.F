      subroutine inf_qnt_c
	  use solid_mod
	  use vertex_mod
      implicit none
      integer :: i,k,ig,indk0,indk,k0,gbs_number
      real*8 :: qq(3)
	  integer, allocatable :: ifla(:)
c -------- Symmetrization k-points into another q - points --------------------
	  allocate(ifla(nqdiv_c))
	  ifla=0
	  i=npnt_c
	  do k=1,npnt_c
	    i_kref_c(k)=k
	    k_group_c(k)=1
	    call zone1_number(pnt_c(:,k),rb0_c,ndiv_c,indk0)
	    index_k1_c(indk0)=k
	    ifla(indk0)=1
	    do ig=1,ngroup
          call rotate(pnt_c(1,k),pnt_c(2,k),pnt_c(3,k),qq(1),qq(2),
     &                qq(3),u(2,ig),2)   !! direct operation: qq=Rk
	      if(ifia(ig).eq.-1) qq=-qq  !! ???????????????
	      call zone1_number(qq,rb0_c,ndiv_c,indk)
	      if(indk>0.and.ifla(indk).eq.0) then
	        i=i+1
	        i_kref_c(i)=k
	        k_group_c(i)=ig
	        index_k1_c(indk)=i
	        pnt_c(:,i)=qq
	        ifla(indk)=1
	      endif
        enddo   !! over igroup
	  enddo   !! over k
	  deallocate(ifla)
	  do ig=1,ngroup
	    do k=1,nqdiv_c
          call rotate(pnt_c(1,k),pnt_c(2,k),pnt_c(3,k),qq(1),qq(2),
     &                qq(3),u(2,ig),2)
	      call zone1_number(qq,rb0_c,ndiv_c,indk)
	      if(indk>0) then
	        do i=1,nqdiv_c
	          if(index_k1_c(indk)==i) then
	            ip_k_c(i,ig)=k
	            exit
	          endif
	        enddo
	      endif
	    enddo
	  enddo
	  k_star_c=0
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
	    k_star_c(k0)=k_star_c(k0)+1
	  enddo
	  k_star_max_c=maxval(k_star_c)
	  allocate(k_list_c(k_star_max_c,npnt_c))
	  k_star_c=0
	  do k=1,nqdiv_c
	    k0=i_kref_c(k)
	    k_star_c(k0)=k_star_c(k0)+1
	    k_list_c(k_star_c(k0),k0)=k
	  enddo
c --------------- K_SYM_0_C --------------------------------
	  allocate(k_sym_0_c(nqdiv_c,ngroup))
	  allocate(g_sym_0_c(nqdiv_c,ngroup))
      do k0=1,nqdiv_c
	    do ig=1,ngroup
	      call rotate(pnt_c(1,k0),pnt_c(2,k0),pnt_c(3,k0),qq(1),qq(2),
     &                qq(3),u(2,ig),2)  !! direct operation = A*k0
	      call zone1_number(qq,rb0_c,ndiv_c,k)
	      k=index_k1_c(k)
	      k_sym_0_c(k,ig)=k0
	      qq=pnt_c(:,k)-qq
          g_sym_0_c(k,ig)=gbs_number(qq)
	    enddo
	  enddo
      end
