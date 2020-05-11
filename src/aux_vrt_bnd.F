      subroutine aux_vrt_bnd(z_red,a_red,phf,ax)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, intent(out) :: z_red(nfun_red,nbndf_bnd,nqdiv_c,
     &                                 nspin),
     &                           a_red(nbasmpw_red,nbndf_bnd,nqdiv_c,
     &                                 nspin),
     &                           phf(nr_full_red,nqdiv_c),
     &                           ax(nrel*nplw_gw_red,nbndf_bnd,nqdiv_c,
     &                              nspin)
      integer :: k,ispin,i,i1,isort,k0,ir,ind,ka,ka0,irl,itr,km1,km1f,
     &           lm1,km,kmf,lm,ii,iii,indx,icg1,icg2,icg,ie,km2,km0,
     &           lget,l2,lm2,id,istg,istgg,j,jnd,n,nb,nn,gbs_number
      real*8 :: v(3),phase,pi2
      real*8, allocatable :: fifim(:,:,:,:,:)
      complex*16, allocatable :: s(:,:),aa(:,:)
      pi2=pi+pi
c -------- Z_RED, A_RED for K=1,nqdiv_c -------------------------------
      z_red=(0.d0,0.d0)
      a_red=(0.d0,0.d0)
	  do ispin=1,nspin
	    do k=1,npnt_c
          ka=k_a_from_c(k)
	      z_red(:,:,k,ispin)=z_bnd_red(:,:,ka,ispin)
	      a_red(:,:,k,ispin)=a_bnd_red(:,:,ka,ispin)
	    enddo
        do k=npnt_c+1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(z_red(1,1,k,ispin),z_red(1,1,k0,ispin),
     &               k_group_c(k),pnt_c(1,k),nfun_red,
     &               n_low_bnd(ka0,ispin),maxb_red,io_lem_red,lmb_red,
     &               indbasa_red,limlb_red,1)
	      call sym_a(a_red(1,1,k,ispin),a_red(1,1,k0,ispin),ka,
     &               k_group_c(k),pnt_c(1,k),nbasmpw_red,
     &               n_low_bnd(ka0,ispin),nbask_red(ka0),
     &               indgb_red(1,ka),iplf_bk_red(1,ka0))
        enddo
      enddo
      if(irel/=2.or..not.rel_interst) 
     &	call check_band_states_red(z_red,a_red)
c --------------------------------------------------------------
      do k=1,nqdiv_c
        do ir=1,nr_full_red
          v=r0b0_red(:,1)*r_pnt_red(1,ir)+r0b0_red(:,2)*r_pnt_red(2,ir)
     &                                   +r0b0_red(:,3)*r_pnt_red(3,ir)
          phase=pi2*dot_product(pnt_c(:,k),v)
	      phf(ir,k)=dcmplx(cos(phase),sin(phase))
        enddo
      enddo
      allocate(nind_fif(maxel_red,maxel_red,nsort))
      nindm_fif=0
      do itr=1,2
        do isort=1,nsort
          do km1=1,lfunm_red(isort)
            km1f=lme_red(km1,isort)
            lm1=lm_isz(km1f,isort)
            do km=1,lfunm_red(isort)
              kmf=lme_red(km,isort)
              lm=lm_isz(kmf,isort)
              if(irel<=1) then
                ii = max0(lm,lm1)
                iii = min0(lm,lm1)
                indx = (ii*(ii-1))/2 + iii
                icg1 = indxcg(indx)
                icg2 = indxcg(indx+1) - 1 
              else
                icg1=indxcgr(lm,lm1,1,0)
                icg2=indxcgr(lm,lm1,2,0)
              endif
              ind=0
              do icg = icg1, icg2
                if(irel<=1) lm2 = jcg(icg)
                if(irel==2) lm2 = jcgr(icg)
                l2=lget(lm2)
                if(l2>lmpb_red(isort)) cycle
                do ie=1,ntle_pb_red(l2,isort)
                  km2=ind_prod_red(ie,lm2,isort)
                  ind=ind+1
                  if(itr==2) ind_fif(ind,km,km1,isort)=km2
                enddo
              enddo
              if(itr==1) then
                nind_fif(km,km1,isort)=ind
                nindm_fif=max(nindm_fif,ind)
              endif
            enddo
          enddo
        enddo   !! over isort
        if(itr==1) allocate(ind_fif(nindm_fif,maxel_red,maxel_red,
     &                      nsort))
      enddo   !! over itr
      allocate(fif0(nrel,nindm_fif,maxel_red,maxel_red,nsort,nspin))
      fif0=0.d0
c ---------------------------------------------------------------            
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do isort=1,nsort
          do i1=1,lfunm_red(isort)
            do i=1,lfunm_red(isort)
              do km=1,nind_fif(i,i1,isort)
                km0=ind_fif(km,i,i1,isort)
                do irl=1,nrel
                  fif0(irl,km,i,i1,isort,ispin)=
     &                             fifim(irl,i,i1,km0,isort)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
      ax=(0.d0,0.d0)
      allocate(s(nplw_gw_red,nbasmpw_red/nrel))
      allocate(aa(nbasmpw_red/nrel,nbndf_bnd))
      do k=1,nqdiv_c
        ka=k_a_from_c(k)
        k0=i_kref_c(k)
        ka0=k_a_from_c(k0)
        n=nbask_red(ka0)/nrel
        nn=nplwgw_red(ka0)
        s=(0.d0,0.d0)
        do i=1,nn
          ind=indpw_gw_red(i,ka)
          do j=1,n
            jnd=indgb_red(j,ka)
            v=gbs(:,jnd)-gbs(:,ind)
            id=gbs_number(v)
            if(complex_ro) then
              s(i,j)=dcmplx(sovr(id),sovi(id))
            else
              s(i,j)=dcmplx(sovr(id),0.d0)
            endif
          enddo
        enddo
        do ispin=1,nspin
          nb=n_low_bnd(ka0,ispin)
          do irl=1,nrel
            istg=(irl-1)*n
            istgg=(irl-1)*nn+1
            aa(1:n,1:nb)=a_red(istg+1:istg+n,1:nb,k,ispin)
            call zgemm('n','n',nn,nb,n,(1.d0,0.d0),s,nplw_gw_red,
     &                 aa,nbasmpw_red/nrel,(0.d0,0.d0),
     &                 ax(istgg,1,k,ispin),nrel*nplw_gw_red)
          enddo
        enddo
      enddo
      deallocate(aa,s)
      call timel('**** AUX_VRT_BND finished **********')
      end
      
      
      
      
      subroutine aux_vrt_loc(fif,z_red)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      real*8, intent(out) :: fif(nrel,maxel_red**2,n_pbmtm_red,nsort,
     &                           nspin)
      complex*16, intent(out) :: z_red(nfun_red,nbndf_bnd,nqdiv_c,
     &                                 nspin)
      integer :: k,ispin,i,i1,isort,ij,k0,kf,ka,ka0,irl
      real*8 :: pi2
      real*8, allocatable :: fifim(:,:,:,:,:)
      pi2=pi+pi
c ---------------------------------------------------------------            
      allocate(fifim(nrel,maxel_red,maxel_red,n_pbmtm_red,nsort))
      do ispin=1,nspin
        call fi0_full_red(fifim,ispin,ispin)
        do isort=1,nsort
          ij=0
          do i1=1,lfunm_red(isort)
            do i=1,lfunm_red(isort)
              ij=ij+1
              do irl=1,nrel
                fif(irl,ij,:,isort,ispin)=fifim(irl,i,i1,:,isort)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(fifim)
c -------- Z_RED for K=1,nqdiv_c -------------------------------
      z_red=(0.d0,0.d0)
	  do ispin=1,nspin
	    do k=1,npnt_c
          kf=k_a_from_c(k)
	      z_red(:,:,k,ispin)=z_bnd_red(:,:,kf,ispin)
	    enddo
        do k=npnt_c+1,nqdiv_c
	      k0=i_kref_c(k)
          ka=k_a_from_c(k)
          ka0=k_a_from_c(k0)
          call sym_z(z_red(1,1,k,ispin),z_red(1,1,k0,ispin),
     &               k_group_c(k),pnt_c(1,k),nfun_red,
     &               n_low_bnd(ka0,ispin),maxb_red,io_lem_red,lmb_red,
     &               indbasa_red,limlb_red,1)
        enddo
      enddo
      call timel('**** AUX_VRT_LOC finished **********')
      end
