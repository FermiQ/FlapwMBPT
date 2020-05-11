      subroutine tsig_rs_to_ks_mm(t_rs,t_pw,k1,z_red,ispin)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: k1,ispin
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          t_rs(n_pbmt_red,nfun_red,nqdiv_c)
      complex*16, intent(inout) :: t_pw(n_pbtot_red,nbndf_bnd,nqdiv_c)
      integer :: nk,kq,kq0,nkq,kfq0,iq,iq0,iqa0,nk0
      real*8 :: v(3)
      complex*16, allocatable :: tm(:,:),tt(:,:)
      allocate(tm(nfun_red,nbndf_bnd))
	  do iq=1,nqdiv_c
	    iq0=i_kref_c(iq)
        iqa0=k_a_from_c(iq0)
	    nk=nplwgw_red(iqa0)
        nk0=n_pbmt_red+nk
	    v=pnt_c(:,k1)-pnt_c(:,iq)
	    call zone1_number(v,rb0_c,ndiv_c,kq)
	    kq=index_k1_c(kq)     !! for K'-Q
	    kq0=i_kref_c(kq)
        kfq0=k_a_from_c(kq0)
	    nkq=n_low_bnd(kfq0,ispin)
        allocate(tt(n_pbmt_red,nkq))
	    tm(:,1:nkq)=conjg(z_red(:,1:nkq,kq))
        call zgemm('n','n',n_pbmt_red,nkq,nfun_red,(1.d0,0.d0),
     &	           t_rs(1,1,iq),n_pbmt_red,tm,nfun_red,(0.d0,0.d0),
     &             tt,n_pbmt_red)
        t_pw(1:n_pbmt_red,1:nkq,iq)=t_pw(1:n_pbmt_red,1:nkq,iq)+tt
        deallocate(tt)
      enddo   !! over iq
      deallocate(tm)
      end
