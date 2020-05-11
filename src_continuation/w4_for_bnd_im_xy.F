      subroutine w4_for_bnd_im_xy(vim,wim,w_tau,v2_q,fif,pw_bos)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	real*8, intent(in) :: fif(maxel_red**2,n_pbmtm_red,nsort,nspin),
     &	                  w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,
     &                            nqdiv_c)
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c),
     &                          pw_bos(nr_full_red,nplw_gw_red,nqdiv_c)
	real*8, intent(out) :: vim(nr_full_red,maxel_red,maxel_red,natom,
     &	                       nqdiv_c),
     &                       wim(nr_full_red,n_pbmt_red,nqdiv_c,
     &                           0:n_tau/2)
	integer :: iatom,isort,ind_tau,ir,iq0,ij,j,i,n1,i_tau,mx2,
     &           n11,ind1,iq,k11,ii
      real*8, allocatable :: tmp22(:,:),tmp33(:,:)
	complex*16, allocatable :: v2r(:,:,:),tmp(:,:),tt(:,:)
c ----------------------------------------------------------------
      wim=0.d0
      vim=0.d0
      mx2=maxel_red**2
      allocate(v2r(nqdiv_c,nr_full_red,n_pbmt_red))
      allocate(tmp22(nr_full_red,n_pbmtm_red))
      allocate(tmp33(nr_full_red,mx2))
      allocate(tt(n_pbtot_red,n_pbtot_red))
      allocate(tmp(nr_full_red,n_pbmt_red))
c ---- Static part ------------------------------------------------	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
	  call zgemm('n','n',nr_full_red,n_pbmt_red,nplwgw_red(iq0),
     &	         (1.d0,0.d0),pw_bos(1,1,iq),nr_full_red,
     &             v2_q(n_pbmt_red+1,1,iq),n_pbtot_red,(0.d0,0.d0),tmp,
     &             nr_full_red)
	  call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	  v2r(k11,:,:)=tmp
      enddo   !! over iq
	call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &          v2r,1)
	v2r=v2r/nqdiv_c
	do iatom=1,natom
	  isort=is(iatom)
	  n1=lfunm_red(isort)
	  n11=n_pbmt0_red(isort)
	  ind1=iopb_red(iatom)-1
        do ir=1,nqdiv_c
	    do j=1,nr_full_red
	      do i=1,n11
	        tmp22(j,i)=v2r(ir,j,ind1+i)
	      enddo
	    enddo
	    call dgemm('n','t',nr_full_red,n1*n1,n11,1.d0,
     &               tmp22,nr_full_red,fif(1,1,isort,nspin),mx2,
     &               0.d0,tmp33,nr_full_red)
          ij=0
          do j=1,n1
            do i=1,n1
              ij=ij+1
              vim(:,i,j,iatom,ir)=tmp33(:,ij)
            enddo
          enddo
        enddo   !! over ir
      enddo  !! over iatom
c ---- Dynamic part ------------------------------------------------
      do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
c --------------------------------------------------------------- 
        v2r=(0.d0,0.d0)
        do iq=1,nqdiv_c
          call boson_unpack_tau(tt,n_pbtot_red,
     &                          w_tau(1,1,ind_tau,iq),n_pbtot_red,
     &                          n_pbtot_red)
          iq0=i_kref_c(iq)
	    call zgemm('n','n',nr_full_red,n_pbmt_red,nplwgw_red(iq0),
     &	           (1.d0,0.d0),pw_bos(1,1,iq),nr_full_red,
     &               tt(n_pbmt_red+1,1),n_pbtot_red,(0.d0,0.d0),tmp,
     &               nr_full_red)
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	    v2r(k11,:,:)=tmp
        enddo   !! over iq
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &            v2r,1)
	  v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
	    do j=1,n_pbmt_red
	      do i=1,nr_full_red
	        wim(i,j,ir,i_tau)=v2r(ir,i,j)
	      enddo
	    enddo
        enddo   !! over ir
      enddo  !! over ind_tau
      deallocate(tmp22,tmp33,tt,v2r,tmp)
	if(nproc_t/=1) then
	  ii=n_pbmt_red*nr_full_red*nqdiv_c*(n_tau/2+1)
	  call dgop(wim,ii,'  +',comm_t)
	endif
      end