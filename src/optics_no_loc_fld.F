      subroutine optics_no_loc_fld(nrax,wrx,e)
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
      complex*16, intent(out) :: e(3,3,0:nrax)
	  integer :: k,ispin,n,ib,i,j,ibl,ind_k,jb,jbl,ig,i_nu,i_len
      real*8 :: wg,pintra(3,3),pr1(3,3),pr2(3,3),prr(3,3),de,
     &          fermi_dirac,fk,f,fk1,f1,pi8
      complex*16 :: p1(3,3),p2(3,3),per(3,3),cc
      real*8, allocatable :: loss(:,:)
      complex*16, allocatable :: ea(:,:,:),p_inter(:,:,:),ew(:),
     &                           e_intra(:,:,:),e_inter(:,:,:)
      pi8=8.d0*pi
      allocate(ew(0:nrax))
      allocate(p_inter(3,3,0:nrax))
      pintra=0.d0
      p_inter=(0.d0,0.d0)
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
                do j=1,3
                  do i=1,3
                    pintra(i,j)=pintra(i,j)+f*ea(ib,ib,i)*ea(ib,ib,j)
                  enddo
                enddo
              else if(ibl/=jbl) then
                de=e_bnd(jb,k,ispin)-chem_pot
                fk1=fermi_dirac(de)
                de=e_bnd(jb,k,ispin)-e_bnd(ib,k,ispin)
                f=wg*(fk1-fk)/de
                do i_nu=0,nrax
                  ew(i_nu)=de*de-wrx(i_nu)**2
                enddo
                do j=1,3
                  do i=1,3
                    f1=f*ea(ib,jb,i)*conjg(ea(ib,jb,j))
                    do i_nu=0,nrax
                      p_inter(i,j,i_nu)=p_inter(i,j,i_nu)+f1/ew(i_nu)
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
	  if(nproc_k/=1) then
        call DGOP(p_inter,18*(nrax+1),'  +',comm_k)
        call DGOP(pintra,9,'  +',comm_k)
      endif
      if(nspin==1) then
        pintra=2.d0*pintra
        p_inter=2.d0*p_inter
      endif
c ---------------- Symmetrization ------------------------------
      prr=0.d0
      do ig=1,ngroup
        do j=1,3
          call rotate_r(pintra(j,:),pr1(j,:),u(2,ig),1)
        enddo
        do j=1,3
          call rotate_r(pr1(1,j),pr2(1,j),u(2,ig),1)
        enddo
        prr=prr+pr2/ngroup
      enddo
      pintra=prr
      do i=0,nrax
        per=(0.d0,0.d0)
        do ig=1,ngroup
          do j=1,3
            call rotate_c(p_inter(j,:,i),p1(j,:),u(2,ig),1)
          enddo
          do j=1,3
            call rotate_c(p1(1,j),p2(1,j),u(2,ig),1)
          enddo
          per=per+p2/ngroup
        enddo
        p_inter(:,:,i)=per
      enddo
c --------------------------------------------------------------
      allocate(e_intra(3,3,0:nrax))
      allocate(e_inter(3,3,0:nrax))
      e_intra=(0.d0,0.d0)
      e_inter=(0.d0,0.d0)
      do i_nu=0,nrax
        cc=pi8/wrx(i_nu)**2
        do j=1,3
          do i=1,3
            if(i==j) e_intra(i,i,i_nu)=(1.d0,0.d0)
            e_intra(i,j,i_nu)=e_intra(i,j,i_nu)-cc*pintra(i,j)
            e_inter(i,j,i_nu)=e_inter(i,j,i_nu)-pi8*p_inter(i,j,i_nu)
          enddo
        enddo
      enddo
      deallocate(p_inter)
      e=e_intra+e_inter
c ------------ Loss function -------------------------------------
      allocate(loss(0:nrax,3))
      do i_nu=0,nrax
        loss(i_nu,1)=-imag((1.d0,0.d0)/e_intra(1,1,i_nu))
        loss(i_nu,2)=-imag((1.d0,0.d0)/e_inter(1,1,i_nu))
        loss(i_nu,3)=-imag((1.d0,0.d0)/
     &                 (e_intra(1,1,i_nu)+e_inter(1,1,i_nu)))
      enddo
c ----------------------------------------------------------------
	  i_len=len_trim(allfile)
	  if(maswrk) then
        write(2,*)' Plasmon frequency = ',sqrt(pi8*pintra(1,1))*evolt/2,
     &            ' eV'
	    open(3,file=allfile(1:i_len)//'_opt_eps_Re')
        write(3,*)'             INTRA_BAND   INTER_BAND     FULL'
	    do i_nu=0,nrax
          f=real(e_intra(1,1,i_nu))
          f1=real(e_inter(1,1,i_nu))
 	      write(3,'(f11.5,3(1x,e12.5))')real(wrx(i_nu))*evolt/2,f,f1,
     &                                  f+f1
	    enddo
        close(3)
	    open(3,file=allfile(1:i_len)//'_opt_eps_Im')
        write(3,*)'             INTRA_BAND   INTER_BAND     FULL'
	    do i_nu=0,nrax
          f=imag(e_intra(1,1,i_nu))
          f1=imag(e_inter(1,1,i_nu))
 	      write(3,'(f11.5,3(1x,e12.5))')real(wrx(i_nu))*evolt/2,f,f1,
     &                                  f+f1
	    enddo
        close(3)
	    open(3,file=allfile(1:i_len)//'_opt_loss')
        write(3,*)'             INTRA_BAND   INTER_BAND     FULL'
	    do i_nu=0,nrax
 	      write(3,'(f11.5,3(1x,e12.5))')real(wrx(i_nu))*evolt/2,
     &       loss(i_nu,1),loss(i_nu,2),loss(i_nu,3)
	    enddo
        close(3)
      endif
      deallocate(e_intra,e_inter,loss)
      end
