      subroutine k_sym_intact
c     Finds symmetry operations which left k-point intact 
	  use manager_mod
	  use parallel_mod
	  use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: k,ig,ii,gbs_number,i_len,i,iv(3),k0,ka
      real*8 :: v(3),s,pi2a
      pi2a=(pi+pi)/par
	  allocate(num_sym_k(nqdiv))
	  allocate(list_sym_k(ngroup,nqdiv))
	  allocate(gtild_sym_k(ngroup,nqdiv))
	  allocate(n_class_k(npnt))
	  allocate(ind_class_k(ngroup,npnt))
	  do k=1,nqdiv
	    num_sym_k(k)=0
	    do ig=1,ngroup
	      if(ip_k(k,ig)==k) then
	        num_sym_k(k)=num_sym_k(k)+1
	        list_sym_k(num_sym_k(k),k)=ig
            call rotate(pnt(1,k),pnt(2,k),pnt(3,k),v(1),v(2),v(3),
     &                  u(2,ig),2)
            gtild_sym_k(num_sym_k(k),k)=gbs_number(pnt(:,k)-v)
	      endif
	    enddo
        if(k<=npnt) call classes_k(k)
	  enddo
	  allocate(num_sym_k0(npnt))
	  allocate(list_sym_k0(ngroup,npnt))
	  do k=1,npnt
	    num_sym_k0(k)=0
	    do ii=1,num_sym_k(k)
	      ig=list_sym_k(ii,k)
          call rotate(pnt(1,k),pnt(2,k),pnt(3,k),v(1),v(2),v(3),u(2,ig),
     &                2)
          s=(pnt(1,k)-v(1))**2+(pnt(2,k)-v(2))**2+(pnt(3,k)-v(3))**2
	      if(s<1.d-8) then
	        num_sym_k0(k)=num_sym_k0(k)+1
	        list_sym_k0(num_sym_k0(k),k)=ig
	      endif
	    enddo
	  enddo
      if(maswrk) then
	    i_len=len_trim(allfile)
	    open(3,file=allfile(1:i_len)//'_K_points')
	    write(3,'(a42)')' Reciprocal space translations (in 2pi/A):'
        write(3,'(a5,3(1x,e15.8),a3,f12.7)')' Q_A ',(gbas(i,1),i=1,3),
     &    ' = ',sqrt(dot_product(gbas(:,1),gbas(:,1)))
        write(3,'(a5,3(1x,e15.8),a3,f12.7)')' Q_B ',(gbas(i,2),i=1,3),
     &    ' = ',sqrt(dot_product(gbas(:,2),gbas(:,2)))
        write(3,'(a5,3(1x,e15.8),a3,f12.7)')' Q_C ',(gbas(i,3),i=1,3),
     &    ' = ',sqrt(dot_product(gbas(:,3),gbas(:,3)))
	    write(3,'(a10)')' K-points:'
	    do k=1,nqdiv
	      k0=i_kref(k)
	      call zone1_number_v(pnt(:,k),rb0,ndiv,iv)
          s=sqrt(dot_product(pnt(:,k),pnt(:,k)))*pi2a
	      write(3,'(i4,3(1x,e15.8),3i4,1x,a1,1x,e12.5,a5)')k,
     &     (pnt(i,k),i=1,3),(iv(i),i=1,3),sym_label(k0),s,' a.u.'
        enddo
        close(3)
	    open(3,file=allfile(1:i_len)//'_K_points_C')
	    write(3,'(a42)')' Reciprocal space translations (in 2pi/A):'
	    write(3,'(a5,3(1x,e15.8))')' Q_A ',(gbas(i,1),i=1,3)
	    write(3,'(a5,3(1x,e15.8))')' Q_B ',(gbas(i,2),i=1,3)
	    write(3,'(a5,3(1x,e15.8))')' Q_C ',(gbas(i,3),i=1,3)
	    write(3,'(a10)')' K-points:'
	    do k=1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k0)
	      call zone1_number_v(pnt_c(:,k),rb0_c,ndiv_c,iv)
          s=sqrt(dot_product(pnt_c(:,k),pnt_c(:,k)))*pi2a
	      write(3,'(i4,3(1x,e15.8),3i4,1x,a1,1x,e12.5,a5)')k,
     &     (pnt_c(i,k),i=1,3),(iv(i),i=1,3),sym_label(ka),s,' a.u.'
        enddo
        close(3)
      endif
      end
