      subroutine w4_for_bnd_mi_xy(vmi,wmi,w_tau,v2_q,fif,pw_bos)
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
	real*8, intent(out) :: vmi(maxel_red,maxel_red,natom,nr_full_red,
     &	                       nqdiv_c),
     &                       wmi(n_pbmt_red,nr_full_red,nqdiv_c,
     &                           0:n_tau/2)
	integer :: iatom,isort,ind_tau,ir,iq0,ij,j,i,n1,i_tau,mx2,
     &           n11,ind1,iq,k11,ii
      real*8, allocatable :: tmp22(:,:),tmp33(:,:)
	complex*16, allocatable :: v2r(:,:,:),tmp(:,:),tt(:,:)
c ----------------------------------------------------------------
      wmi=0.d0
      vmi=0.d0
      mx2=maxel_red**2
      allocate(v2r(nqdiv_c,n_pbmt_red,nr_full_red))
      allocate(tmp22(n_pbmtm_red,nr_full_red))
      allocate(tmp33(mx2,nr_full_red))
      allocate(tt(n_pbtot_red,n_pbtot_red))
      allocate(tmp(n_pbmt_red,nr_full_red))
c ---- Static part ------------------------------------------------	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
        iq0=i_kref_c(iq)
	  call zgemm('n','c',n_pbmt_red,nr_full_red,nplwgw_red(iq0),
     &	         (1.d0,0.d0),v2_q(1,n_pbmt_red+1,iq),n_pbtot_red,
     &             pw_bos(1,1,iq),nr_full_red,(0.d0,0.d0),tmp,
     &             n_pbmt_red)
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
	        tmp22(i,j)=v2r(ir,ind1+i,j)
	      enddo
	    enddo
	    call dgemm('n','n',n1*n1,nr_full_red,n11,1.d0,
     &	           fif(1,1,isort,1),mx2,tmp22,n_pbmtm_red,0.d0,tmp33,
     &               mx2)
          ij=0
          do j=1,n1
            do i=1,n1
              ij=ij+1
              vmi(i,j,iatom,:,ir)=tmp33(ij,:)
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
	    call zgemm('n','c',n_pbmt_red,nr_full_red,nplwgw_red(iq0),
     &	           (1.d0,0.d0),tt(1,n_pbmt_red+1),n_pbtot_red,
     &               pw_bos(1,1,iq),nr_full_red,(0.d0,0.d0),tmp,
     &               n_pbmt_red)
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	    v2r(k11,:,:)=tmp
        enddo   !! over iq
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red*nr_full_red,
     &            v2r,1)
	  v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
	    do j=1,nr_full_red
	      do i=1,n_pbmt_red
	        wmi(i,j,ir,i_tau)=v2r(ir,i,j)
	      enddo
	    enddo
        enddo   !! over ir
      enddo  !! over ind_tau
      deallocate(tmp22,tmp33,tt,v2r,tmp)
	if(nproc_t/=1) then
	  ii=n_pbmt_red*nr_full_red*nqdiv_c*(n_tau/2+1)
	  call dgop(wmi,ii,'  +',comm_t)
	endif
      end