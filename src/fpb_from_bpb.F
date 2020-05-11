      subroutine fpb_from_bpb(pl,pf,k)
	use atom_mod
	use parallel_mod
	use solid_mod
	use vertex_mod
	implicit none
	integer, intent(in) :: k
	complex*16, intent(in) :: pl(npb_g2_w2_max,npb_g2_w2_max)
	complex*16, intent(out) :: pf(n_pbtot,n_pbtot)
	integer :: n_pbt,np,k1
	complex*16, allocatable :: tmp(:,:)
	pf=(0.d0,0.d0)
      k1=i_kref(k)
	n_pbt=n_pbmt+nplwgw(k1)
	np=npb_g2_w2(k1)
	allocate(tmp(n_pbt,np))
	call zgemm('n','n',n_pbt,np,np,(1.d0,0.d0),
     &	         m_p_overl_s(1,1,k),n_pbtot,pl,npb_g2_w2_max,
     &           (0.d0,0.d0),tmp,n_pbt)
	call zgemm('n','c',n_pbt,n_pbt,np,(1.d0,0.d0),tmp,n_pbt,
     &	         m_p_overl_s(1,1,k),n_pbtot,(0.d0,0.d0),pf,n_pbtot)
      deallocate(tmp)
      end