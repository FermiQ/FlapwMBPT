      subroutine opt_chi0_q0(b11,b12,b21,b22,nne,nrax,npb,vie,wrx)
	  use atom_mod
	  use manager_mod
      use models_mod
	  use parallel_mod
	  use units_mod
	  use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nne,nrax,npb
      complex*16, intent(in) :: vie(npb,nne),wrx(0:nrax)
      complex*16, intent(out) :: b11(0:nrax,3,3),b21(nne,0:nrax,3),
     &                           b12(nne,0:nrax,3),b22(nne,nne,0:nrax)
	  integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,k0,i_star,inu,ipb
      real*8 :: de,fermi_dirac,fk,fk1,pi8,gtild(3,2)
      complex*16 :: cc,ew,c1
      real*8, allocatable :: wfk(:),fifim(:,:,:,:)
      complex*16, allocatable :: ea(:,:,:),zn(:,:),
     &                           an(:,:),ppm(:,:,:),ppe(:,:,:),
     &                           bb11(:,:),bb12(:,:),bb21(:,:),
     &                           bb22(:,:),pp2(:,:,:)
      pi8=8.d0*pi
      b11=(0.d0,0.d0)
      b12=(0.d0,0.d0)
      b21=(0.d0,0.d0)
      b22=(0.d0,0.d0)
      gtild=0.d0
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
      allocate(bb11(3,3))
      allocate(bb12(nne,3))
      allocate(bb21(nne,3))
      allocate(bb22(nne,nne))
      do ispin=1,nspin
        call fi0_full(fifim,ispin,ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k0=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k0,ispin)
          allocate(zn(nfun,n))
          allocate(an(nbasmpw,n))
          allocate(ea(n,n,3))
          allocate(wfk(n))
          allocate(ppm(n,n,npb))
          allocate(ppe(n,n,nne))
          allocate(pp2(n,n,nne))
          do ib=1,n
            de=e_bnd(ib,k0,ispin)-chem_pot
            fk=fermi_dirac(de)
            wfk(ib)=betta_t*fk*(1.d0-fk)
          enddo
          do i_star=1,k_star(k0)
            k=k_list(i_star,k0)
            call sym_z(zn,z_bnd(1,1,ind_k,ispin),k_group(k),pnt(1,k),
     &                 nfun,n,maxb,io_lem,lmb,indbasa,limlb,0)
            call sym_a(an,ev_bnd(1,1,ind_k,ispin),k,k_group(k),pnt(1,k),
     &                 nbasmpw,n,nbask(k0),indgb(1,k),iplf_bk(1,k0))
            call bands_gradient(n,n,k,ispin,zn,an,ea)
            ea=(0.d0,-2.d0)*ea
            call ppm_factors(k,k,1,fifim,zn,zn,an,an,ppm,gtild,npb,n,n,
     &                       0,nbndf)
            call zgemm('n','n',n*n,nne,npb,(1.d0,0.d0),ppm,n*n,vie,npb,
     &                 (0.d0,0.d0),ppe,n*n)
c ----------------------------------------------------------------
            do inu=0,nrax
              bb11=(0.d0,0.d0)
              bb12=(0.d0,0.d0)
              bb21=(0.d0,0.d0)
              bb22=(0.d0,0.d0)
              pp2=(0.d0,0.d0)
              do ib=1,n
                ibl=ind_block_k(ib,k0,ispin)
                de=e_bnd(ib,k0,ispin)-chem_pot
                fk=fermi_dirac(de)
                do jb=1,n
                  jbl=ind_block_k(jb,k0,ispin)
                  if(ib==jb) then
                    do j=1,3
                      do i=1,3
                        bb11(i,j)=bb11(i,j)+wfk(ib)*ea(ib,ib,i)
     &                                     *ea(ib,ib,j)/wrx(inu)**2
                      enddo
                    enddo
                  else if(ibl/=jbl) then
                    de=e_bnd(jb,k0,ispin)-chem_pot
                    fk1=fermi_dirac(de)
                    de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                    ew=de*de-wrx(inu)**2
                    ew=(fk1-fk)/ew
                    do j=1,3
                      do i=1,3
                        c1=conjg(ea(ib,jb,i))*ea(ib,jb,j)
                        bb11(i,j)=bb11(i,j)+c1*ew/de
                      enddo
                    enddo
                    do ipb=1,nne
                      cc=ew*ppe(ib,jb,ipb)
                      do i=1,3
                        bb12(ipb,i)=bb12(ipb,i)-cc*conjg(ea(ib,jb,i))
                      enddo
                      cc=ew*conjg(ppe(ib,jb,ipb))
                      do i=1,3
                        bb21(ipb,i)=bb21(ipb,i)-cc*ea(ib,jb,i)
                      enddo
                    enddo
                    cc=ew*de
                    do j=1,nne
                      pp2(ib,jb,j)=cc*ppe(ib,jb,j)
                    enddo
                  endif
                enddo   !! over jb
              enddo    !! over ib
              call zgemm('c','n',nne,nne,n*n,(1.d0,0.d0),ppe,n*n,pp2,
     &                   n*n,(0.d0,0.d0),bb22,nne)
              do j=1,3
                do i=1,3
                  b11(inu,i,j)=b11(inu,i,j)+bb11(i,j)/nqdiv/amega
                enddo
                b12(:,inu,j)=b12(:,inu,j)+bb12(:,j)/nqdiv/sqrt(amega)
                b21(:,inu,j)=b21(:,inu,j)+bb21(:,j)/nqdiv/sqrt(amega)
              enddo
              b22(:,:,inu)=b22(:,:,inu)+bb22/nqdiv
            enddo    !! over inu
          enddo   !! over i_star
          deallocate(ea,zn,an,wfk,ppm,ppe,pp2)
        enddo   !! over ind_k
      enddo    !! over ispin
      deallocate(fifim,bb11,bb12,bb21,bb22)
	  if(nproc_k/=1) then
        call DGOP(b11,18*(nrax+1),'  +',comm_k)
        call DGOP(b12,6*nne*(nrax+1),'  +',comm_k)
        call DGOP(b21,6*nne*(nrax+1),'  +',comm_k)
        call DGOP(b22,2*nne*nne*(nrax+1),'  +',comm_k)
      endif
      if(nspin==1) then
        b11=2.d0*b11
        b12=2.d0*b12
        b21=2.d0*b21
        b22=2.d0*b22
      endif
      end
