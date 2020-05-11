      subroutine p0_head_tensor(nrax,wrx,pab)
	use atom_mod
	use manager_mod
      use models_mod
	use parallel_mod
	use units_mod
	use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nrax
      complex*16, intent(in) :: wrx(0:nrax)
      complex*16, intent(out) :: pab(3,3,0:nrax)
	integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,ig,inu
      real*8 :: wg,de,fermi_dirac,fk,f,fk1,f1,pi8
      complex*16 :: p1(3,3),p2(3,3),per(3,3),cc
      complex*16, allocatable :: ea(:,:,:),ew(:)
      pi8=8.d0*pi
      allocate(ew(0:nrax))
      pab=(0.d0,0.d0)
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          allocate(ea(n,n,3))
          call bands_gradient(nbndf,n,k,ispin,z_bnd(1,1,ind_k,ispin),
     &                        ev_bnd(1,1,ind_k,ispin),ea)
          ea=(0.d0,-2.d0)*ea
          do ib=1,n
            ibl=ind_block_k(ib,k,ispin)
            de=e_bnd(ib,k,ispin)-chem_pot
            fk=fermi_dirac(de)
            wg=wgt(k)/amega
            do jb=1,n
              jbl=ind_block_k(jb,k,ispin)
              if(ib==jb) then
                f=betta_t*wg*fk*(1.d0-fk)
                do inu=0,nrax
                  cc=wrx(inu)**2
                  do j=1,3
                    do i=1,3
                      pab(i,j,inu)=pab(i,j,inu)
     &                            +f*ea(ib,ib,i)*ea(ib,ib,j)/cc
                    enddo
                  enddo
                enddo
              else if(ibl/=jbl) then
                de=e_bnd(jb,k,ispin)-chem_pot
                fk1=fermi_dirac(de)
                de=e_bnd(jb,k,ispin)-e_bnd(ib,k,ispin)
                f=wg*(fk1-fk)/de
                do inu=0,nrax
                  ew(inu)=de*de-wrx(inu)**2
                enddo
                do j=1,3
                  do i=1,3
                    f1=f*conjg(ea(ib,jb,i))*ea(ib,jb,j)
                    do inu=0,nrax
                      pab(i,j,inu)=pab(i,j,inu)+f1/ew(inu)
                    enddo
                  enddo
                enddo
              endif
            enddo
          enddo
          deallocate(ea)
        enddo
      enddo
      deallocate(ew)
	if(nproc_k/=1) call DGOP(pab,18*(nrax+1),'  +',comm_k)
      if(nspin==1) pab=2.d0*pab
c ---------------- Symmetrization ------------------------------
      do i=0,nrax
        per=(0.d0,0.d0)
        do ig=1,ngroup
          do j=1,3
            call rotate_c(pab(j,:,i),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_c(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        pab(:,:,i)=per
      enddo
      end