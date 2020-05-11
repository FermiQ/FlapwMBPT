      subroutine prod_bas_red
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer :: limax,ndim0,in0,jn0,in,jn,isort,l2,li,li1,nl0,l,k,ie0,
     &           le0,mtpb0,ie,le,mtpb,i,ind_k,n_pbt,kb,ir,ia,ib,ic,ii,
     &           ind,j,jj
      real*8 :: dqdall
	  integer, allocatable :: nl(:,:),ind_l(:,:,:)
      real*8, allocatable :: work(:),v_tmpr(:,:)
	  complex*16, allocatable :: v_tmp(:,:),s(:,:),b(:,:)
	  limax=nrel*maxb_red+1
c	  key=0 - Run to find the dimensions
	  allocate(nl(limax,limax))
	  allocate(ind_l(limax,limax,maxpb+1))
	  call prod_b1(ndim0,nl,ind_l,limax)
	  if(.not.allocated(lmpb_red)) allocate(lmpb_red(nsort))
	  lmpb_red=-1
	  do isort=1,nsort
	    do l2=0,lmpb(isort)
          if(lim_pb_mt_red(l2,isort)==0) cycle
	      do in0=1,lfun_red(isort)
	        in=le_red(in0,isort)
 		    li=li_val(in,isort)
	        do jn0=1,lfun_red(isort)
	          jn=le_red(jn0,isort)
			  li1=li_val(jn,isort)
	          do nl0=1,nl(li,li1)
	            if(ind_l(li,li1,nl0).eq.l2)
     &	          lmpb_red(isort)=max(lmpb_red(isort),l2)
	          enddo
	        enddo   !! over jn0
	      enddo   !! over in0
	    enddo   !! over l2
	  enddo   !! over isort
	  maxpb_red=maxval(lmpb_red)
      if(maswrk) write(iun,*)' Max dimension in reduced PB search is ',
     &                       ndim0
	  if(.not.allocated(lfun_pb_red)) allocate(lfun_pb_red(nsort))
	  if(.not.allocated(n_pbmt0_red)) allocate(n_pbmt0_red(nsort))
	  if(.not.allocated(iopb_red)) allocate(iopb_red(natom))
	  if(.not.allocated(iopb0_red)) allocate(iopb0_red(nsort))
	  if(.not.allocated(ntle_pb_red))
     &	allocate(ntle_pb_red(0:maxpb_red,nsort))
      ntle_pb_red=0
	  call prod_b2_red(ndim0,0,nl,ind_l,limax)
	  if(allocated(ind_pb_red)) deallocate(ind_pb_red)
	  allocate(ind_pb_red(maxlfpb_red,nsort))
	  if(allocated(lval_pb_red)) deallocate(lval_pb_red)
	  allocate(lval_pb_red(maxlfpb_red,nsort))
	  if(allocated(fun_pb_red)) deallocate(fun_pb_red)
	  allocate(fun_pb_red(maxpbr_red))
	  if(allocated(dfun_pb_red)) deallocate(dfun_pb_red)
	  allocate(dfun_pb_red(maxpbr_red))
	  if(allocated(lm_pbmt_red)) deallocate(lm_pbmt_red)
	  allocate(lm_pbmt_red(n_pbmtm_red,nsort))
	  if(allocated(lf_pb_red)) deallocate(lf_pb_red)
	  allocate(lf_pb_red(n_pbmtm_red,nsort))
	  if(allocated(ind_prod0_red)) deallocate(ind_prod0_red)
	  allocate(ind_prod0_red(max_e_pb_red,0:maxpb,nsort))
	  if(allocated(ind_prod_red)) deallocate(ind_prod_red)
	  allocate(ind_prod_red(max_e_pb_red,limlpb,nsort))
c	  key=1 - Run to fill the arrays
	  call prod_b2_red(ndim0,1,nl,ind_l,limax)
	  if(allocated(fi0_red)) deallocate(fi0_red)
      allocate(fi0_red(maxlfun_red,maxlfun_red,maxlfpb_red,nsort,nspin,
     &                 nspin))
	  call get_fi0_red
	  if(allocated(fi0_f_red)) deallocate(fi0_f_red)
      allocate(fi0_f_red(maxlfun,maxlfun,maxlfpb_red,nsort,nspin,nspin))
	  call get_fi0_f_red
	  deallocate(nl,ind_l)
	  n_pbtot_red=n_pbmt_red+nplw_gw_red
	  nrs_red=nfun_red+nrel*nr_full_red
	  nra_red=nfun_red+nbasmpw_red
	  npb_rs_red=max(n_pbmt_red,nr_full_red)
	  np_rs_red=n_pbmt_red+nr_full_red
	  do k=1,npnt_c
        kb=k_a_from_c(k)
	    i=n_pbmt_red+nplwgw_red(kb)
        call size_shift_par(i,nproc_b,ndim4_pbr(1,k),n4_mpi_pbr(1,k))
	  enddo
      call size_shift_par(n_pbmt_red,nproc_b,ndim4_pbr_mt,
     &                    n4_mpi_pbr_mt)
	  if(maswrk) write(iun,'(a30,i6,a5,i6,a7,i6)')
     &	             'Reduced Product Basis: MT-part',n_pbmt_red,
     &                 ' Int:',nplw_gw_red,' Total:',n_pbtot_red
      if(.not.allocated(ndim3_nrs_red2)) 
     &                        allocate(ndim3_nrs_red2(nproc_k))
      if(.not.allocated(n3_mpi_nrs_red2)) 
     &                        allocate(n3_mpi_nrs_red2(nproc_k))
      if(.not.allocated(ndim_nrs_red2)) 
     &                        allocate(ndim_nrs_red2(nproc))
      if(.not.allocated(n_mpi_nrs_red2)) 
     &                        allocate(n_mpi_nrs_red2(nproc))
      call size_shift_par(nrs_red**2,nproc,ndim_nrs_red2,n_mpi_nrs_red2)
	  call size_shift_par(nrs_red**2,nproc_k,ndim3_nrs_red2,
     &	                  n3_mpi_nrs_red2)
c ---------- Overlaps with the full Product Basis --------------------	
	  if(allocated(rpb_fpb)) deallocate(rpb_fpb)
	  allocate(rpb_fpb(max_e_pb_red,max_e_pb,0:maxpb_red,nsort))
	  rpb_fpb=0.d0
	  allocate(work(0:maxnrad))
	  do isort=1,nsort
	    do l=0,lmpb_red(isort)
	      do ie0=1,ntle_pb_red(l,isort)
	        le0=ind_prod0_red(ie0,l,isort)
	        mtpb0=ind_pb_red(le0,isort)
	        do ie=1,ntle_pb(l,isort)
	          le=ind_prod0(ie,l,isort)
	          mtpb=ind_pb(le,isort)
	          do ir=0,nrad(isort)
                work(ir)=fun_pb_red(ir+mtpb0)*fun_pb(ir+mtpb)
     &	  	            *dr(ir,isort)*r(ir,isort)**2
	          enddo
              rpb_fpb(ie0,ie,l,isort)=dqdall(h(isort),work,nrad(isort))
	        enddo   !! over ie
	      enddo    !! over ie0
	    enddo   !! over l
	  enddo  !! over isort
	  deallocate(work)
      if(allocated(rpb_fpb_int)) deallocate(rpb_fpb_int)
      allocate(rpb_fpb_int(nplw_gw,nplw_gw_red,npnt))
      allocate(s(nplw_gw_red,nplw_gw_red))
      allocate(b(nplw_gw,nplw_gw_red))
      rpb_fpb_int=(0.d0,0.d0)
      do k=1,npnt
        if(nplwgw_red(k)/=0) call s_overlap(1,s,nplw_gw_red,
     &                                      nplwgw_red(k),
     &                                      indpw_gw_red(1,k))
        s=s*amega
        do j=1,nplwgw_red(k)
          jj=indpw_gw_red(j,k)
          do i=1,nplwgw(k)
            ii=indpw_gw(i,k)
            ia=igbs(1,jj)-igbs(1,ii)
            ib=igbs(2,jj)-igbs(2,ii)
            ic=igbs(3,jj)-igbs(3,ii)
            ind=indplw(ia,ib,ic)
            if(complex_ro) then
              b(i,j)=dcmplx(sovr(ind),sovi(ind))
            else
              b(i,j)=dcmplx(sovr(ind),0.d0)
            endif
          enddo
        enddo
        if(nplwgw_red(k)/=0) then
          call zgemm('n','n',nplwgw(k),nplwgw_red(k),nplwgw_red(k),
     &               (1.d0,0.d0),b,nplw_gw,s,nplw_gw_red,(0.d0,0.d0),
     &               rpb_fpb_int(1,1,k),nplw_gw)
        endif
      enddo   !! over ind_k
      deallocate(s,b)
	  call trans_pb_pw_red
c ------------------------------------------------------------------
	  if(mesh_k_c) then
        if(.not.allocated(v_red_q))
     &	  allocate(v_red_q(n_pbtot_red,n_pbtot_red,npnt))
        if(.not.allocated(w_nu_q_b))
     &	  allocate(w_nu_q_b(n_pbtot_red,n_pbtot_red,ndim3_nu,
     &                        ndim3_k(me_k+1)))
        v_red_q=(0.d0,0.d0)
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
	      n_pbt=n_pbmt+nplwgw(k)
	      allocate(v_tmp(n_pbt,n_pbt))
          if(ncmpl==1) then
            allocate(v_tmpr(n_pbt,n_pbt))
            call v_coul_full_r(nplwgw(k),v_tmpr,ind_k,0)
            call pb_c_from_r(n_pbt,n_pbt,n_pbt,v_tmp,v_tmpr,pnt(1,k))
            deallocate(v_tmpr)
          else
            call v_coul_full(nplwgw(k),v_tmp,ind_k,0)
          endif
          call rpb_from_fpb_full(v_red_q(1,1,k),v_tmp,k,n_pbtot_red,
     &                           n_pbt)
          deallocate(v_tmp)
        enddo
        if(nproc_k/=1) call dgop(v_red_q,2*n_pbtot_red**2*npnt,'  +',
     &                           comm_k)
        call coulomb_ev_red
      endif
      call bz_interp_prep_red
      call timel('**** prod_bas_red finished *********')
      end
