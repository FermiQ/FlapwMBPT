      subroutine sym_v_mpb(iq,n0,w)
	use atom_mod
	use manager_mod
	use parallel_mod
	use units_mod
	use solid_mod
      implicit none
	integer, intent(in) :: n0,iq
	complex*16, intent(inout)  :: w(n0,n0)
	integer :: ig,j,i,ii,j0,j1,i0,i1,gbs_number,igtild,n
	real*8 :: pi2,v(3),phase
	complex*16 :: cc
	complex*16, allocatable :: tmp(:,:),tmi(:,:),w1(:,:)
	n=nplwgw(iq)
	pi2=pi+pi
      allocate(w1(n0,n0))
	allocate(tmp(n_pbmt,n_pbmt))
	allocate(tmi(n_pbmt,nplw_gw))
      w1=(0.d0,0.d0)
	do ii=1,num_sym_k(iq)
	  ig=list_sym_k(ii,iq)
	  igtild=g_sym_0(iq,ig)
        tmp=w(1:n_pbmt,1:n_pbmt)
c ----------------- MM part ---------------------------------		
	  call sym_w_mt(ig,pnt(1,iq),tmp,n_pbmt)
        w1(1:n_pbmt,1:n_pbmt)=w1(1:n_pbmt,1:n_pbmt)+tmp
c ----------------- MI part ---------------------------------
        tmi(:,1:n)=w(1:n_pbmt,n_pbmt+1:n_pbmt+n)
        call sym_w_mi(iq,ig,tmi)
        w1(1:n_pbmt,n_pbmt+1:n_pbmt+n)=w1(1:n_pbmt,n_pbmt+1:n_pbmt+n)
     &                                +tmi(:,1:n)
c ----------------- II part ---------------------------------
        do j=1,n
	    j0=indpw_gw(j,iq)  !! G' in GBS-list
	    v=gbs(:,j0)+gbs(:,igtild)  !! G'+G_A
          j1=gbs_number(v)   !! G'+G_A in GBS-list
          j1=gbs_sym(j1,ig)   !! A^(-1)*(G'+G_A) in GBS-list
          j1=iplf_gk(j1,iq)
          do i=1,n
	      i0=indpw_gw(i,iq)
	      v=gbs(:,i0)+gbs(:,igtild)  !! G+G_A
            i1=gbs_number(v)   !! G+G_A in GBS-list
            i1=gbs_sym(i1,ig)   !! A^(-1)*(G+G_A) in GBS-list
            i1=iplf_gk(i1,iq)
	      v=gbs(:,j0)-gbs(:,i0)
	      phase=pi2*dot_product(v,shift(:,ig))
	      cc=dcmplx(cos(phase),sin(phase))
	      cc=cc*w(n_pbmt+i1,n_pbmt+j1)
            w1(n_pbmt+i,n_pbmt+j)=w1(n_pbmt+i,n_pbmt+j)+cc
	    enddo
	  enddo   !! over j
      enddo
      w1=w1/num_sym_k(iq)
      do j=1,n
        do i=1,n_pbmt
          w1(j,i)=conjg(w1(i,j))
        enddo
      enddo
      w=w1
	deallocate(tmp,tmi,w1)
      end