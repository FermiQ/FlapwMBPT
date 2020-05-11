      subroutine output_pdos_x
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: iedos,ispin,i_len,isort,li,ind_dos,l,ie
      real*8 :: deltae,e0,consp(2),om0,fct
      integer, allocatable :: nd_dos(:),nmpi_dos(:)
      real*8, allocatable :: weigd(:,:,:),tm(:),
     &	                   weigh(:,:,:),pdos(:,:,:,:),pdos_at(:,:,:),
     &                       pdos_diag(:,:,:,:,:)
      consp(1)=-1.d0
      consp(nspin)=1.d0
      fct=evolt/2
      allocate(weigd(nbndf,npnt,nspin),weigh(nbndf,npnt,nspin))
      allocate(pdos(nrel*maxb+1,nsort,nspin,0:ndos))
      allocate(pdos_at(0:nsort,nspin,0:ndos))
      allocate(pdos_diag(maxntle,nrel*maxb+1,nsort,nspin,0:ndos))
      allocate(tm(nbndf))
      deltae=(emax_pdos_x-emin_pdos_x)/ndos
      allocate(nd_dos(nproc_t))
      allocate(nmpi_dos(nproc_t))
      call size_shift_par(ndos+1,nproc_t,nd_dos,nmpi_dos)
      pdos=0.d0
      pdos_at=0.d0
      pdos_diag=0.d0
      do ind_dos=1,nd_dos(me_t+1)
        iedos=nmpi_dos(me_t+1)+ind_dos-1
        om0=emin_pdos_x+iedos*deltae
        e0=om0+chem_pot
        call sum_up_weights_x(weigd,om0,nbndf)	    
        call pdos_x(pdos(1,1,1,iedos),pdos_at(0,1,iedos),
     &              pdos_diag(1,1,1,1,iedos),weigd,nbndf)
      enddo
      if(nproc_t/=1) then
        call dgop(pdos,(nrel*maxb+1)*nsort*nspin*(1+ndos),'  +',
     &	        comm_t)
        call dgop(pdos_at,(nsort+1)*nspin*(1+ndos),'  +',
     &	        comm_t)
        call dgop(pdos_diag,maxntle*(nrel*maxb+1)*nsort*nspin*(1+ndos),
     &            '  +',comm_t)
      endif
      pdos=pdos/fct
      pdos_at=pdos_at/fct
      pdos_diag=pdos_diag/fct
      deallocate(nd_dos,nmpi_dos)
      if(maswrk) then
        i_len=len_trim(allfile)
c ---- Partial DOSs --------------------------------------------
        do isort=1,nsort    
          if(ubi=='dft') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.dft')	      
          if(ubi==' hf') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.hft')	      
          if(ubi==' gw') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.gw')	      
          if(ubi==' qp') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.qp')	      
          if(ubi=='psi') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.psi')	      
          if(ubi=='bsp') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x.bsp')
          do iedos=0,ndos
            om0=emin_pdos_x+iedos*deltae
            write(3,'(9f8.3)')om0*fct,
     &			  ((pdos(li,isort,ispin,iedos)*consp(ispin)
     &                                          *nhsort(isort),
     &                  li=1,nrel*min(3,lmb(isort))+1),ispin=1,nspin)
          enddo
          close(3)
        enddo   !! over isort
c ---- Partial DOS_DIAGs --------------------------------------------
        do isort=1,nsort      
          if(ubi=='dft') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.dft')	      
          if(ubi==' hf') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.hft')	      
          if(ubi==' gw') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.gw')
          if(ubi==' qp') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.qp')
          if(ubi=='psi') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.psi')
          if(ubi=='bsp') open(3,file=allfile(1:i_len)
     &	                    //txtel(isort)(1:2)//'_pdos_x_diag.bsp')
          do li=1,nrel*lmb(isort)+1
            write(3,*)'----------- Li = ',li
            if(irel.ne.2) then
              l=li-1
            else if(irel.eq.2) then
              l=li/2
            endif
            do iedos=0,ndos
              om0=emin_pdos_x+iedos*deltae
              write(3,'(9f8.3)')om0*fct,
     &			  ((pdos_diag(ie,li,isort,ispin,iedos)*consp(ispin)
     &                                          *nhsort(isort),
     &                  ie=1,min(4,ntle(l,isort))),ispin=1,nspin)
            enddo
          enddo
          close(3)
        enddo   !! over isort
c ---- Partial DOS_ATs --------------------------------------------
        if(ubi=='dft') open(3,file=allfile(1:i_len)//'_pdos_x_at.dft')
        if(ubi==' hf') open(3,file=allfile(1:i_len)//'_pdos_x_at.hft')
        if(ubi==' gw') open(3,file=allfile(1:i_len)//'_pdos_x_at.gw')
        if(ubi==' qp') open(3,file=allfile(1:i_len)//'_pdos_x_at.qp')
        if(ubi=='psi') open(3,file=allfile(1:i_len)//'_pdos_x_at.psi')
        if(ubi=='bsp') open(3,file=allfile(1:i_len)//'_pdos_x_at.bsp')
        do iedos=0,ndos
          om0=emin_pdos_x+iedos*deltae
          write(3,'(11f7.2)')om0*fct,
     &			  (pdos_at(0,ispin,iedos)*consp(ispin),
     &               (pdos_at(isort,ispin,iedos)*consp(ispin)
     &                                          *nhsort(isort),
     &                  isort=1,min(4,nsort)),ispin=1,nspin)
        enddo
        close(3)
      endif
      deallocate(weigd,weigh,pdos,tm,pdos_at,pdos_diag)
      call timel('***** OUTPUT_PDOS_X finished *******')
      end