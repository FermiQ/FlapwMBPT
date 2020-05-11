      subroutine p0_q0(b11,b12,b21,b22,nne,npb,vie,wp)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use units_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: nne,npb
      complex*16, intent(in) :: vie(npb,nne)
      real*8, intent(out) :: wp(3,3)
      complex*16, intent(out) :: b11(3,3,ndim3_nu),b21(nne,3,ndim3_nu),
     &                           b12(nne,3,ndim3_nu),
     &                           b22(nne,nne,ndim3_nu)
      integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,k0,i_star,ipb,ind_nu,
     &           i_nu
      real*8 :: de,fermi_dirac,fk,fk1,pi8,gtild(3,2),wp0(3,3),wn
      complex*16 :: cc,ew,c1
      real*8, allocatable :: wfk(:),fifim(:,:,:,:)
      complex*16, allocatable :: ea(:,:,:),zn(:,:),
     &                           an(:,:),ppm(:,:,:),ppe(:,:,:),
     &                           pp2(:,:,:)
      pi8=8.d0*pi
      b11=(0.d0,0.d0)
      b12=(0.d0,0.d0)
      b21=(0.d0,0.d0)
      b22=(0.d0,0.d0)
      wp=0.d0
      gtild=0.d0
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
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
            wp0=0.d0
            do ib=1,n
              de=e_bnd(ib,k0,ispin)-chem_pot
              fk=fermi_dirac(de)
              do j=1,3
                do i=1,3
                  wp0(i,j)=wp0(i,j)-wfk(ib)*ea(ib,ib,i)
     &                                     *ea(ib,ib,j)
                enddo
              enddo
            enddo    !! over ib
            do j=1,3
              do i=1,3
                wp(i,j)=wp(i,j)+wp0(i,j)/nqdiv/amega
              enddo
            enddo
            do ind_nu=1,ndim3_nu
              i_nu=me_t*ndim3_nu+ind_nu-1
              wn=w_nu(i_nu)
              pp2=(0.d0,0.d0)
              do ib=1,n
                ibl=ind_block_k(ib,k0,ispin)
                de=e_bnd(ib,k0,ispin)-chem_pot
                fk=fermi_dirac(de)
                do jb=1,n
                  jbl=ind_block_k(jb,k0,ispin)
                  if(ibl/=jbl) then
                    de=e_bnd(jb,k0,ispin)-chem_pot
                    fk1=fermi_dirac(de)
                    de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                    ew=de*de+wn*wn
                    ew=(fk1-fk)/ew
                    do j=1,3
                      do i=1,3
                        c1=conjg(ea(ib,jb,i))*ea(ib,jb,j)
                        b11(i,j,ind_nu)=b11(i,j,ind_nu)+c1*ew/de
                      enddo
                    enddo
                    do ipb=1,nne
                      cc=ew*ppe(ib,jb,ipb)
                      do i=1,3
                        b12(ipb,i,ind_nu)=b12(ipb,i,ind_nu)
     &                                   -cc*conjg(ea(ib,jb,i))
                      enddo
                      cc=ew*conjg(ppe(ib,jb,ipb))
                      do i=1,3
                        b21(ipb,i,ind_nu)=b21(ipb,i,ind_nu)
     &                                   -cc*ea(ib,jb,i)
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
     &                   n*n,(1.d0,0.d0),b22(1,1,ind_nu),nne)
            enddo   !! over ind_nu
          enddo   !! over i_star
          deallocate(ea,zn,an,wfk,ppm,ppe,pp2)
        enddo   !! over ind_k
      enddo    !! over ispin
      deallocate(fifim)
      b11=b11/nqdiv/amega
      b12=b12/nqdiv/sqrt(amega)
      b21=b21/nqdiv/sqrt(amega)
      b22=b22/nqdiv
      if(nproc_k/=1) then
        call DGOP(wp,9,'  +',comm_k)
        call DGOP(b11,18*ndim3_nu,'  +',comm_k)
        call DGOP(b12,6*nne*ndim3_nu,'  +',comm_k)
        call DGOP(b21,6*nne*ndim3_nu,'  +',comm_k)
        call DGOP(b22,2*nne*nne*ndim3_nu,'  +',comm_k)
      endif
      if(nspin==1) then
        wp=2.d0*wp
        b11=2.d0*b11
        b12=2.d0*b12
        b21=2.d0*b21
        b22=2.d0*b22
      endif
      end
