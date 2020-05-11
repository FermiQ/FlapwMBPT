      subroutine chi0_lda_static_1q(iq,isp1,isp2,chi0,n)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: isp1,isp2,n,iq
      complex*16, intent(out) :: chi0(n,n)
      integer :: k,ind_k,i,j,jj,ig,n1,n2,npr,me1,n12,j0,i_star,
     &           iqs,kq,kq0,ii
      real*8 :: gtild(3,2),de,green0,aa,v(3)
      integer, allocatable :: n2d(:),n20(:)
      real*8, allocatable :: gb(:,:),fifim(:,:,:,:)
      complex*16, allocatable :: ps3(:,:,:),tmp(:,:),tmp1(:,:),t(:,:,:),
     &                           zb0(:,:,:),ev0(:,:,:),zb(:,:),ev(:,:)
      npr=nproc_t*nproc_b
      me1=me_t_b+1
      allocate(n2d(npr))
      allocate(n20(npr))
      allocate(tmp(n,n))            
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
      call fi0_full(fifim,isp1,isp2)
      allocate(gb(nbndf,2))
      chi0=(0.d0,0.d0)
      gtild(:,1)=0.d0
      allocate(zb0(nfun,nbndf,npnt))
      allocate(ev0(nbasmpw,nbndf,npnt))
      zb0=(0.d0,0.d0)
      ev0=(0.d0,0.d0)
	  do ind_k=1,ndim3_k(me_k+1)
	    k=n3_mpi_k(me_k+1)+ind_k
        zb0(:,:,k)=z_bnd(:,:,ind_k,isp2)
        ev0(:,:,k)=ev_bnd(:,:,ind_k,isp2)
      enddo
      if(nproc_k/=1) then
        call dgop(zb0,2*nfun*nbndf*npnt,'  +',comm_k)
        call dgop(ev0,2*nbasmpw*nbndf*npnt,'  +',comm_k)
      endif
      allocate(zb(nfun,nbndf))
      allocate(ev(nbasmpw,nbndf))
c ---------------------------------------------------------------      
      do i_star=1,k_star(iq)
        iqs=k_list(i_star,iq)
        tmp=(0.d0,0.d0)
	    do ind_k=1,ndim3_k(me_k+1)
	      k=n3_mpi_k(me_k+1)+ind_k
          v=pnt(:,k)-pnt(:,iqs)
          call zone1_number(v,rb0,ndiv,kq)
          kq=index_k1(kq)    !! for K-Q
          kq0=i_kref(kq)
          gtild(:,2)=v-pnt(:,kq)
          n1=n_bnd(k,isp1)
          n2=n_bnd(kq0,isp2)
          do i=1,n1
            de=e_bnd(i,k,isp1)-chem_pot
            gb(i,1)=green0(de,betta_t)
          enddo
          do i=1,n2
            de=e_bnd(i,kq0,isp2)-chem_pot
            gb(i,2)=green0(de,betta_t)
          enddo
          call size_shift_par(n2,npr,n2d,n20)
          allocate(ps3(n1,n2d(me1),n))
          call sym_z(zb,zb0(1,1,kq0),k_group(kq),pnt(1,kq),nfun,n2,maxb,
     &               io_lem,lmb,indbasa,limlb,0)
	      call sym_a(ev,ev0(1,1,kq0),kq,k_group(kq),pnt(1,kq),nbasmpw,
     &               n2,nbask(kq0),indgb(1,kq),iplf_bk(1,kq0))
          call ppm_factors0(k,kq,iqs,fifim,z_bnd(1,1,ind_k,isp1),zb,
     &                      ev_bnd(1,1,ind_k,isp1),ev,ps3,gtild,n,n1,
     &                      n2d(me1),n20(me1),nbndf)
          allocate(t(n1,n2d(me1),n))
          do j=1,n2d(me1)
            j0=n20(me1)+j  
            do i=1,n1
              de=e_bnd(j0,kq0,isp2)-e_bnd(i,k,isp1)
              if(abs(de)<1.d-8) then
                aa=betta_t*(1.d0+gb(i,1))*gb(i,1)
              else
                aa=(gb(i,1)-gb(j0,2))/de
              endif
              do jj=1,n
                t(i,j,jj)=aa*ps3(i,j,jj)
              enddo
            enddo
          enddo
          t=wgt(k)*t
          n12=n1*n2d(me1)
          call zgemm('c','n',n,n,n12,(1.d0,0.d0),ps3,n12,t,n12,
     &               (1.d0,0.d0),tmp,n)
          deallocate(ps3,t)
        enddo   !! over ind_k
        if(nproc_k/=1) call dgop(tmp,2*n*n,'  +',comm_k)
        if(npr/=1) call dgop(tmp,2*n*n,'  +',comm_t_b)
c --------- Symmetrization --------------------------------------
        allocate(tmp1(n,n))
  	    do ii=1,num_sym_k(iqs)
	      ig=list_sym_k(ii,iqs)
          tmp1=tmp
          call sym_chi_q(ig,tmp1,n,iqs)
          chi0=chi0+tmp1
        enddo
        deallocate(tmp1)
      enddo   !! iver i_star
      chi0=chi0/ngroup
      deallocate(fifim,gb,n2d,n20,tmp,zb0,ev0,zb,ev)
      call timel('**** CHI0_LDA_STATIC_1Q finished ***')
      end
