      subroutine bands_vertex1
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,k,ibnd,ib,ind_k,ii,ii0,iatom,ind,isort,i,j,
     & 	       j0,lf,ie,lm,li,l,km1,nbas,nbasr,irl,iii,jj
      if(allocated(z_bnd_red)) deallocate(z_bnd_red)
      allocate(z_bnd_red(nfun_red,nbndf_bnd,npnt,nspin))
      if(allocated(a_bnd_red)) deallocate(a_bnd_red)
      allocate(a_bnd_red(nbasmpw_red,nbndf_bnd,npnt,nspin))
      z_bnd_red=(0.d0,0.d0)
      a_bnd_red=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          nbasr=nbask_red(k)/nrel
          nbas=nbask(k)/nrel
          do ibnd=1,n_low_bnd(k,ispin)
            ib=ind_bands_bnd(ibnd,k,ispin)
            ii=0
            do iatom=1,natom
              isort=is(iatom)
              ii0=io_lem(iatom)-1
              do i=1,lfunm_red(isort)
                ii=ii+1
                ind=lme_red(i,isort) !! == km in full basis
                lf=lf_isz(ind,isort)
                ie=ie_isz(ind,isort)
                lm=lm_isz(ind,isort)
                li=li_val(lf,isort)
                if(irel.le.1) l=li-1
                if(irel.eq.2) l=li/2
                if(correlated(ie,l,isort)=='A') then
                  km1=indbasa(2,ie,lm,isort)
                  z_bnd_red(ii,ibnd,k,ispin)=
     &              z_bnd(ii0+ind,ib,ind_k,ispin)
     &             *augm_coef(1,li,isort,ispin)
     &             +z_bnd(ii0+km1,ib,ind_k,ispin)
     &             *augm_coef(2,li,isort,ispin)
                else if(correlated(ie,l,isort)=='L') then
                  z_bnd_red(ii,ibnd,k,ispin)=
     &              z_bnd(ii0+ind,ib,ind_k,ispin)
                endif
              enddo
            enddo
            do j=1,nbasr  !! in reduced basis for k
              j0=indgb_red(j,k)  !! in GBS-list
              ii=iplf_bk(j0,k)    !! in full basis for k
              do irl=1,nrel
                jj=(irl-1)*nbasr+j
                iii=(irl-1)*nbas+ii
                a_bnd_red(jj,ibnd,k,ispin)=ev_bnd(iii,ib,ind_k,ispin)
              enddo
            enddo
          enddo
        enddo
      enddo
      if(nproc_k/=1) then
        call dgop(a_bnd_red,2*nbasmpw_red*nbndf_bnd*npnt*nspin,'  +',
     &            comm_k)
        call dgop(z_bnd_red,2*nfun_red*nbndf_bnd*npnt*nspin,'  +',
     &            comm_k)
      endif
      if(mesh_k_c) call uu_prepare_ac_old
      end