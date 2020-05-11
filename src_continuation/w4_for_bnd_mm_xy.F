      subroutine w4_for_bnd_mm_xy(v4_r,w4_r,w_tau,v2_q,fif)
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
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
	real*8, intent(out) :: v4_r(maxel_red**2,maxel_red**2,natom**2,
     &	                        nqdiv_c),
     &                       w4_r(n_pbmt_red,n_pbmt_red,nqdiv_c,
     &                            0:n_tau/2)
	integer :: jatom,jsort,iatom,isort,ind_tau,ir,j,i,n1,i_tau,itt,n2,
     &           mx2,n22,jnd1,n11,ind1,jj,i3,i2,ii,i1,jj0,
     &           ii0,iq,k11
      real*8, allocatable :: tmp22(:,:),tmp33(:,:),tt(:,:)
	complex*16, allocatable :: tmp2(:,:),v2r(:,:,:),tmp(:,:)
c ----------------------------------------------------------------
      w4_r=0.d0
      v4_r=0.d0
      mx2=maxel_red**2
      allocate(v2r(nqdiv_c,n_pbmt_red,n_pbmt_red))
      allocate(tmp2(nqdiv_c,n_pbmtm_red**2))
      allocate(tmp22(n_pbmtm_red,n_pbmtm_red))
      allocate(tmp33(mx2,n_pbmtm_red))
      allocate(tt(mx2,mx2))
      allocate(tmp(n_pbmt_red,n_pbmt_red))
c ---- Static part ------------------------------------------------	 
      v2r=(0.d0,0.d0)
      do iq=1,nqdiv_c
	  call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
	  v2r(k11,:,:)=v2_q(1:n_pbmt_red,1:n_pbmt_red,iq)
      enddo   !! over iq
	call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,1)
	v2r=v2r/nqdiv_c
	itt=0
	do jatom=1,natom
	  jsort=is(jatom)
	  n2=lfunm_red(jsort)
	  n22=n_pbmt0_red(jsort)
	  jnd1=iopb_red(jatom)-1
	  do iatom=1,natom
	    itt=itt+1
	    isort=is(iatom)
	    n1=lfunm_red(isort)
	    n11=n_pbmt0_red(isort)
	    ind1=iopb_red(iatom)-1
          do ir=1,nqdiv_c
	      do j=1,n22
	        do i=1,n11
	          tmp22(i,j)=v2r(ir,ind1+i,jnd1+j)
	        enddo
	      enddo
	      call dgemm('n','n',n1*n1,n22,n11,1.d0,fif(1,1,isort,1),mx2,
     &	             tmp22,n_pbmtm_red,0.d0,tmp33,mx2)
	      call dgemm('n','t',n1*n1,n2*n2,n22,1.d0,tmp33,mx2,
     &	             fif(1,1,jsort,nspin),mx2,0.d0,tt,mx2)
            jj=0
            do i3=1,n2
              do i2=1,n1
                jj=jj+1
                ii=0
                do i1=1,n2
                  do i=1,n1
                    ii=ii+1
                    jj0=n2*(i3-1)+i1
                    ii0=n1*(i2-1)+i
                    v4_r(ii,jj,itt,ir)=tt(ii0,jj0)
                  enddo
                enddo
              enddo
            enddo
          enddo   !! over ir
        enddo  !! over iatom
      enddo   !! over jatom
c ---- Dynamic part ------------------------------------------------
      do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
c --------------------------------------------------------------- 
        v2r=(0.d0,0.d0)
        do iq=1,nqdiv_c
	    call zone1_number(pnt_c(1,iq),rb0_c,ndiv_c,k11)
          call boson_unpack_tau(tmp,n_pbmt_red,
     &                          w_tau(1,1,ind_tau,iq),n_pbtot_red,
     &                          n_pbmt_red)
	    v2r(k11,:,:)=tmp
        enddo   !! over iq
	  call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),n_pbmt_red**2,v2r,1)
	  v2r=v2r/nqdiv_c
c --------------------------------------------------------------
        do ir=1,nqdiv_c
	    do j=1,n_pbmt_red
	      do i=1,n_pbmt_red
	        w4_r(i,j,ir,i_tau)=v2r(ir,i,j)
	      enddo
	    enddo
        enddo   !! over ir
      enddo  !! over ind_tau
      deallocate(tmp2,tmp22,tmp33,tt,v2r,tmp)
	if(nproc_t/=1) then
	  ii=n_pbmt_red**2*nqdiv_c*(n_tau/2+1)
	  call dgop(w4_r,ii,'  +',comm_t)
	endif
      end