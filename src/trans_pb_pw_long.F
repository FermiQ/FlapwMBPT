      subroutine trans_pb_pw_long
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
      implicit none
	  integer :: k,i,j,ind,nq0,igg,gbs_number,j0
      real*8 :: v(3),qgg
	  real*8, allocatable :: mt_pw_0(:,:)
      complex*16, allocatable :: s(:,:),b(:),c(:),mt_pw(:)
	  if(allocated(m_pw_long)) deallocate(m_pw_long)
	  allocate(m_pw_long(n_pbtot,n_q_long))
	  if(allocated(m_pw_tild_long)) deallocate(m_pw_tild_long)
	  allocate(m_pw_tild_long(n_pbtot,n_q_long))
	  allocate(s(nplw_gw,nplw_gw))
	  allocate(b(nplw_gw))
	  allocate(c(nplw_gw))
      allocate(mt_pw_0(maxlfpb,nsort))
      allocate(mt_pw(n_pbmt))
      do nq0=1,n_q0_long
        k=q0_long_in_npnt(nq0)
        call s_overlap(1,s,nplw_gw,nplwgw(k),indpw_gw(1,k))
        do i=1,size_q_long(nq0)
          ind=index_q_long(2,i,nq0)  !! in Q_long list
          igg=index_q_long(1,i,nq0)  !! in GBS list
c ------- MT part ---------------------------------------------
          v=pnt(:,k)+gbs(:,igg)
          qgg=sqrt(dot_product(v,v))
          call integral_pbmt_pw_0(mt_pw_0,qgg)
          call integral_pbmt_pw(mt_pw,mt_pw_0,v)
          m_pw_long(1:n_pbmt,ind)=conjg(mt_pw)/sqrt(amega)
          m_pw_tild_long(1:n_pbmt,ind)=conjg(mt_pw)/sqrt(amega)
c ------- Interstitial part -----------------------------------
	      do j=1,nplwgw(k)
	        j0=indpw_gw(j,k)
	        v=gbs(:,j0)-gbs(:,igg)
            j0=gbs_number(v)
	        if(complex_ro) then
		      b(j)=dcmplx(sovr(j0),sovi(j0))
	        else
		      b(j)=dcmplx(sovr(j0),0.d0)
	        endif
          enddo
	      call zgemm('c','n',nplwgw(k),1,nplwgw(k),(1.d0,0.d0),s,
     &               nplw_gw,b,nplw_gw,(0.d0,0.d0),c,nplw_gw)
          m_pw_long(n_pbmt+1:n_pbmt+nplwgw(k),ind)=b(1:nplwgw(k))
          m_pw_tild_long(n_pbmt+1:n_pbmt+nplwgw(k),ind)=c(1:nplwgw(k))
        enddo
	  enddo   !! over ind_k
	  deallocate(s,b,c,mt_pw_0,mt_pw)
      end
