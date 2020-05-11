      subroutine trans_pb_pw_red
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
	  integer :: k,n_pbt,ind_gg,i,ia,ib,ic,ii,j,jj,info,ind
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: s(:,:)
	  if(allocated(pw_pb_red)) deallocate(pw_pb_red)
	  allocate(pw_pb_red(n_pbtot_red,npnt))
      if(allocated(pw_pb_tild_red)) deallocate(pw_pb_tild_red)
      allocate(pw_pb_tild_red(n_pbtot_red,npnt))
      allocate(s(nplw_gw_red,nplw_gw_red),ipiv(nplw_gw_red))
	  ind_gg=1
	  do k=1,npnt
	    n_pbt=n_pbmt_red+nplwgw_red(k)
	    call trans_pb_pw_0_red(k,pw_pb_red(1,k),n_pbt,ind_gg)
        pw_pb_tild_red(1:n_pbt,k)=pw_pb_red(1:n_pbt,k)
        do j=1,nplwgw_red(k)
          jj=indpw_gw_red(j,k)
          do i=1,nplwgw_red(k)
            ii=indpw_gw_red(i,k)
            ia=igbs(1,jj)-igbs(1,ii)
            ib=igbs(2,jj)-igbs(2,ii)
            ic=igbs(3,jj)-igbs(3,ii)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              s(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        if(nplwgw_red(k)/=0) call zgesv(nplwgw_red(k),1,s,nplw_gw_red,
     &                                  ipiv,
     &                                  pw_pb_tild_red(n_pbmt_red+1,k),
     &                                  n_pbtot_red,info)
	  enddo   !! over k
	  deallocate(s,ipiv)
      end
