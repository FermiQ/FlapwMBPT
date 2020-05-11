      subroutine psi_psi_pw_prepare(zw,aw)
      use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, intent(in) :: zw(nfun,nbndf_bnd,nqdiv_c,nspin),
     &                          aw(nbasmpw,nbndf_bnd,nqdiv_c,nspin)
      integer :: k,iatom,isort,ind,jbas,jnd,ibas,j,ndimb,ind0,k10f,
     &           gbs_number,iq,kq,k10,kq10,ispin,jspin,kq10f,kf,kqf
      real*8 :: v(3),gtild(3)
      real*8, allocatable :: fifi_j(:,:,:,:,:,:,:),ql(:,:)
      complex*16, allocatable :: tmp(:,:),s(:,:),tmp1(:,:),ff(:,:,:,:)
      allocate(tmp(maxel,nbndf))
      allocate(tmp1(nbndtm,nbndf))
      allocate(s(nbndtm,nbndtm))
      allocate(ff(maxel,maxel,natom,nc_line))
      allocate(ql(3,nc_line))
      do k=1,nc_line
        j=kline_in_npnt_c(k)
        ql(:,k)=pnt_c(:,j)
      enddo
      psi_psi_pw=(0.d0,0.d0)
      allocate(fifi_j(maxlfun,maxlfun,0:2*maxb,nspin,nspin,nsort,
     &                nc_line))
      call fifi_j_prepare(fifi_j,ql,nc_line)
      do jspin=1,nspin
        do ispin=1,nspin
          call wan_wan_pw_0(fifi_j,ff,ql,nc_line,ispin,jspin)
          do iq=1,nc_line
            do k=1,nqdiv_c
              kf=k_a_from_c(k)
              v=pnt_c(:,k)-ql(:,iq)
	          call zone1_number(v,rb0_c,ndiv_c,kq)  ! 1 BZ
              kq=index_k1_c(kq)    !  0 BZ
              kqf=k_a_from_c(kq)
c ----------------- MT contribution --------------------------------
	          do iatom=1,natom
	            isort=is(iatom)
	            ndimb=lfunm(isort)
	            ind0=io_lem(iatom)
	            call zgemm('n','n',ndimb,n_low_bnd(i_kref(kqf),jspin),
     &	                   ndimb,(1.d0,0.d0),ff(1,1,iatom,iq),maxel,
     &                     zw(ind0,1,kq,jspin),nfun,(0.d0,0.d0),tmp,
     &                     maxel)
	            call zgemm('c','n',n_low_bnd(i_kref(kf),ispin),
     &	                   n_low_bnd(i_kref(kqf),jspin),ndimb,
     &                     (1.d0,0.d0),zw(ind0,1,k,ispin),nfun,tmp,
     &                     maxel,(1.d0,0.d0),
     &                     psi_psi_pw(1,1,k,iq,ispin,jspin),nbndf_bnd)
	          enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
	          gtild=v-pnt_c(:,kq)
c ------ We form the overlap matrix -----------------------------
              k10f=i_kref(kf)
              kq10f=i_kref(kqf)
              k10=i_kref_c(k)
              kq10=i_kref_c(kq)
	          do jbas=1,nbask(kq10f)
                jnd=indgb(jbas,kqf)
                do ibas=1,nbask(k10f)
                  ind=indgb(ibas,kf)
                  v=gbs(:,jnd)-gbs(:,ind)-gtild
                  j=gbs_number(v)
	              if(complex_ro) then
	                s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	              else
	                s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	              endif
	            enddo
	          enddo
              call zgemm('n','n',nbask(k10f),n_low_bnd(kq10f,jspin),
     &                   nbask(kq10f),(1.d0,0.d0),s,nbndtm,
     &                   aw(1,1,kq,jspin),nbasmpw,(0.d0,0.d0),tmp1,
     &                   nbndtm)
	          call zgemm('c','n',n_low_bnd(k10f,ispin),
     &	                 n_low_bnd(kq10f,jspin),
     &	                 nbask(k10f),(1.d0,0.d0),aw(1,1,k,ispin),
     &                   nbasmpw,tmp1,nbndtm,(1.d0,0.d0),
     &                   psi_psi_pw(1,1,k,iq,ispin,jspin),nbndf_bnd)
            enddo   !! over k
          enddo   !! over iq
        enddo   !! over ispin
      enddo   !! over jspin
      deallocate(fifi_j,ff,tmp,s,tmp1,ql)
      psi_psi_pw=psi_psi_pw/sqrt(amega)
      end
      
      
      
      
      subroutine ppw_check(ispin,n,aw,zw,chk,k,k1,q)
      use atom_mod
	  use etot_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: n,k,k1,ispin
      real*8, intent(in) :: q
      complex*16, intent(in) :: zw(nfun,n),aw(nbasmpw,n)
      complex*16, intent(out) :: chk(n,nbndf)
      integer :: iatom,isort,ind,jbas,jnd,ibas,j,ndimb,ind0,k00,
     &           gbs_number,n1,ind_k,ln,k10
      real*8 :: v(3),gtild(3)
      real*8, allocatable :: fifi_j(:,:,:)
      complex*16, allocatable :: tmp(:,:),s(:,:),tmp1(:,:),ff(:,:),
     &                           aa(:,:),zz(:,:)
      k10=i_kref(k1)
      n1=n_bnd(k10,ispin)
      chk=(0.d0,0.d0)
      allocate(aa(nbasmpw,n1))
      allocate(zz(nfun,n1))
      aa=(0.d0,0.d0)
      zz=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me_k+1)
        k00=n3_mpi_k(me_k+1)+ind_k
        if(k00/=k10) cycle
        call sym_z(zz,z_bnd(1,1,ind_k,ispin),k_group(k1),pnt(1,k1),nfun,
     &             n1,maxb,io_lem,lmb,indbasa,limlb,0)
        call sym_a(aa,ev_bnd(1,1,ind_k,ispin),k1,k_group(k1),pnt(1,k1),
     &             nbasmpw,n1,nbask(k10),indgb(1,k1),iplf_bk(1,k10))
      enddo
	  if(nproc_k/=1) then
	    call DGOP(zz,2*nfun*n1,'  +',comm_k)
	    call DGOP(aa,2*nbasmpw*n1,'min',comm_k)
	  endif
c ----------------- MT contribution ----------------------------------
      v=-q
	  do iatom=1,natom
	    isort=is(iatom)
	    ndimb=lfunm(isort)
        allocate(tmp(ndimb,n1))
c ----------------------------------------------------------
        ln=lfun(isort)
        allocate(fifi_j(ln,ln,0:2*lmb(isort)))
        call fifi_j_prep(fifi_j,v,ln,isort,ispin)
        allocate(ff(ndimb,ndimb))
        call wan_wan_pw_3(fifi_j,ff,v,iatom,ndimb,ln)
        deallocate(fifi_j)
c -----------------------------------------------------------------------
	    ind0=io_lem(iatom)
	    call zgemm('n','n',ndimb,n1,ndimb,(1.d0,0.d0),ff,ndimb,
     &             zz(ind0,1),nfun,(0.d0,0.d0),tmp,ndimb)
	    call zgemm('c','n',n,n1,ndimb,(1.d0,0.d0),zw(ind0,1),nfun,tmp,
     &             ndimb,(1.d0,0.d0),chk,n)
        deallocate(tmp,ff)
	  enddo  !! over iatom
c ----------------- Interstitial contribution ------------------------
	  gtild=pnt(:,k)+q-pnt(:,k1)
c ------ We form the overlap matrix -----------------------------
      allocate(s(nbask(k),nbask(k10)))
	  do jbas=1,nbask(k10)
        jnd=indgb(jbas,k1)
        do ibas=1,nbask(k)
          ind=indgb(ibas,k)
          v=gbs(:,jnd)-gbs(:,ind)-gtild
          j=gbs_number(v)
	      if(complex_ro) then
	        s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	      else
	        s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	      endif
	    enddo
	  enddo
      allocate(tmp1(nbask(k),n1))
      call zgemm('n','n',nbask(k),n1,nbask(k10),(1.d0,0.d0),s,nbask(k),
     &           aa,nbasmpw,(0.d0,0.d0),tmp1,nbask(k))
	  call zgemm('c','n',n,n1,nbask(k),(1.d0,0.d0),aw,nbasmpw,tmp1,
     &           nbask(k),(1.d0,0.d0),chk,n)
      deallocate(s,tmp1)
      chk=chk/sqrt(amega)
      end
