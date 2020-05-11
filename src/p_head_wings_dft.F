      subroutine p_head_wings_dft(p_head,p_wing,aa00,aa01,ch)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use units_mod
      use solid_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      real*8, intent(out) :: aa00
      complex*16, intent(out) :: p_head(3,3,ndim3_nu),
     &                           p_wing(3,n_opt_pb,ndim3_nu),
     &                           aa01(n_opt_pb)
      integer :: ispin,n,ib,i,j,ibl,ind_k,jb,jbl,k0,ipb,ig,igi,
     &           ind_nu,i_nu,nne,npb,sg
      real*8 :: de,fermi_dirac,fk,fk1,pi8,gtild(3,2),w2,wg,fkk
      complex*16 :: cc,ew,c1,per(3,3),p1(3,3),p2(3,3)
      real*8, allocatable :: wfk(:),fifim(:,:,:,:)
      complex*16, allocatable :: ea(:,:,:),zn(:,:),pwing(:,:,:),
     &                           an(:,:),ppm(:,:),
     &                           bb11(:,:),bb12(:),aa1(:),aa2(:,:),
     &                           aa3(:,:),aa4(:,:)
      pi8=8.d0*pi
      aa00=0.d0
      p_head=(0.d0,0.d0)
      gtild=0.d0
      allocate(fifim(maxel,maxel,n_pbmtm,nsort))
      npb=n_pbmt+nplwgw(1)
      nne=n_opt_pb
      allocate(pwing(3,npb,ndim3_nu))
      pwing=(0.d0,0.d0)
      allocate(aa1(npb))
      aa1=(0.d0,0.d0)
      allocate(bb11(3,3))
      allocate(bb12(3))
      do ispin=1,nspin
        sg=1
        if(ch=='0Z'.and.ispin==1) sg=-1
        call fi0_full(fifim,ispin,ispin)
        do ind_k=1,ndim3_k(me_k+1)
          k0=n3_mpi_k(me_k+1)+ind_k
          wg=wgt(k0)*sg
          n=n_bnd(k0,ispin)
          allocate(zn(nfun,n))
          allocate(an(nbasmpw,n))
          allocate(ea(n,n,3))
          zn=z_bnd(:,1:n,ind_k,ispin)
          an=ev_bnd(:,1:n,ind_k,ispin)
          if(irel/=2) then
            call bands_gradient(n,n,k0,ispin,zn,an,ea)
            ea=(0.d0,-2.d0)*ea
          else
            call c_alpha_bands(n,k0,zn,an,ea)
          endif
          allocate(wfk(n))
          do ib=1,n
            de=e_bnd(ib,k0,ispin)-chem_pot
            fk=fermi_dirac(de)
            wfk(ib)=betta_t*fk*(1.d0-fk)
          enddo
c ----------------------------- Head -------------------------------
          do ind_nu=1,ndim3_nu
            i_nu=me_t*ndim3_nu+ind_nu-1
            w2=w_nu(i_nu)**2
            bb11=(0.d0,0.d0)
            do ib=1,n
              ibl=ind_block_k(ib,k0,ispin)
              de=e_bnd(ib,k0,ispin)-chem_pot
              fk=fermi_dirac(de)
              do jb=1,n
                jbl=ind_block_k(jb,k0,ispin)
                if(ib==jb) then
                  if(wfk(ib)>1.d-10) then
                    if(i_nu==0) then
                      aa00=aa00-wfk(ib)*wg
                    else
                      do j=1,3
                        do i=1,3
                          bb11(i,j)=bb11(i,j)-wfk(ib)*ea(ib,ib,i)
     &                                       *ea(ib,ib,j)/w2
                        enddo
                      enddo
                    endif
                  endif
                else if(ibl/=jbl) then
                  de=e_bnd(jb,k0,ispin)-chem_pot
                  fk1=fermi_dirac(de)
                  fkk=fk1-fk
                  if(abs(fkk)>1.d-10) then
                    de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                    ew=de*de+w2
                    ew=fkk/ew
                    do j=1,3
                      do i=1,3
                        c1=conjg(ea(ib,jb,i))*ea(ib,jb,j)
                        bb11(i,j)=bb11(i,j)+c1*ew/de
                      enddo
                    enddo
                  endif
                endif
              enddo   !! over jb
            enddo    !! over ib
            do j=1,3
              do i=1,3
                p_head(i,j,ind_nu)=p_head(i,j,ind_nu)+bb11(i,j)*wg
     &                                               /amega
              enddo
            enddo
          enddo    !! over ind_nu
c ----------------------------- Wing -------------------------------
          allocate(ppm(n,n))
          do ipb=1,npb
            if(ipb<=n_pbmt) then
              call ppm_factor_mt(ipb,fifim,zn,zn,ppm,n,n)
            else
              call ppm_factor_int(ipb,k0,k0,1,an,an,ppm,gtild,n,n)
            endif
            do ind_nu=1,ndim3_nu
              i_nu=me_t*ndim3_nu+ind_nu-1
              w2=w_nu(i_nu)**2
              bb12=(0.d0,0.d0)
              do ib=1,n
                ibl=ind_block_k(ib,k0,ispin)
                de=e_bnd(ib,k0,ispin)-chem_pot
                fk=fermi_dirac(de)
                do jb=1,n
                  if(abs(ppm(ib,jb))<1.d-10) cycle
                  jbl=ind_block_k(jb,k0,ispin)
                  if(ib==jb) then
                    if(i_nu==0) then
                      aa1(ipb)=aa1(ipb)-wfk(ib)*ppm(ib,ib)*wg
                    endif
                  else if(ibl/=jbl) then
                    de=e_bnd(jb,k0,ispin)-chem_pot
                    fk1=fermi_dirac(de)-fk
                    if(abs(fk1)>1.d-10) then
                      de=e_bnd(jb,k0,ispin)-e_bnd(ib,k0,ispin)
                      ew=de*de+w2
                      ew=fk1/ew
                      cc=ew*ppm(ib,jb)
                      do i=1,3
                        bb12(i)=bb12(i)-cc*conjg(ea(ib,jb,i))
                      enddo
                    endif
                  endif
                enddo   !! over jb
              enddo    !! over ib
              do j=1,3
                pwing(j,ipb,ind_nu)=pwing(j,ipb,ind_nu)+bb12(j)*wg
     &                                                 /sqrt(amega)
              enddo
            enddo    !! over ind_nu
          enddo  !! over ipb
c ------------------------------------------------------------------
          deallocate(ea,zn,an,wfk,ppm)
        enddo   !! over ind_k
      enddo    !! over ispin
      deallocate(fifim,bb11,bb12)
      aa00=aa00/amega
      aa1=aa1/sqrt(amega)
      if(nproc_k/=1) then
        call DGOP(p_head,18*ndim3_nu,'  +',comm_k)
        call DGOP(pwing,6*npb*ndim3_nu,'  +',comm_k)
        call DGOP(aa00,1,'  +',comm_k)
        call DGOP(aa1,2*npb,'  +',comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(aa00,1,'  +',comm_t)
        call DGOP(aa1,2*npb,'  +',comm_t)
      endif
      if(nspin==1.and.irel/=2) then
        p_head=2.d0*p_head
        pwing=2.d0*pwing
        aa00=2.d0*aa00
        aa1=2.d0*aa1
      endif
c ---------------- Symmetrization ----------------------------------
      allocate(aa2(3,npb))
      allocate(aa3(3,npb))
      allocate(aa4(3,npb))
      do ind_nu=1,ndim3_nu
        per=(0.d0,0.d0)
        aa4=(0.d0,0.d0)
        do ig=1,ngroup
          igi=u_inv(ig)
          do j=1,3
            call rotate_c(p_head(j,:,ind_nu),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_c(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
          do j=1,3
            call sym_m_p_s(aa2(j,:),pwing(j,:,ind_nu),1,1,igi)
          enddo
          do j=1,npb
            call rotate_c(aa2(1,j),aa3(1,j),u(2,ig),1)
          enddo
          aa4=aa4+aa3/ngroup
        enddo
        p_head(:,:,ind_nu)=per
        pwing(:,:,ind_nu)=aa4
      enddo
      aa4=(0.d0,0.d0)
      do ig=1,ngroup
        igi=u_inv(ig)
        call sym_m_p_s(aa2(1,:),aa1,1,1,igi)
        aa4(1,:)=aa4(1,:)+aa2(1,:)/ngroup
      enddo
      aa1=aa4(1,:)
      deallocate(aa2,aa3,aa4)
c --------------------------------------------------------------
      do ind_nu=1,ndim3_nu
        call zgemm('n','c',3,nne,npb,(1.d0,0.d0),pwing(1,1,ind_nu),3,
     &             bm_coef,nne,(0.d0,0.d0),p_wing(1,1,ind_nu),3)
      enddo
      call zgemm('n','c',1,nne,npb,(1.d0,0.d0),aa1,1,bm_coef,nne,
     &           (0.d0,0.d0),aa01,1)
      deallocate(pwing,aa1)
      end
