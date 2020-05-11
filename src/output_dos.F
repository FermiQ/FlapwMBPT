      subroutine output_dos
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2 :: s1
      character*3 :: s2
      character*5 :: ss1
      character*6 :: sss
      integer :: iedos,ispin,i_len,isort,li,k,i,ind0,l,ndim,ind_dos,ie,
     &           ind_k,n
      real*8 :: deltae,e0,consp(2),om0,fct,ds1,ds2
      integer, allocatable :: nd_dos(:),nmpi_dos(:)
      real*8, allocatable :: weigd(:),dossp(:,:),pdos(:,:,:,:),
     &                       porb(:,:,:,:),dos_k(:,:,:)
      complex*16, allocatable :: achain(:,:,:)
      allocate(dossp(nspin,0:ndos))
      allocate(dos_k(0:ndos,npnt,nspin))
      allocate(pdos(nrel*maxb+1,nsort,nspin,0:ndos))
      if(irel/=2.and.iprn/=0) allocate(porb(maxel_red,nsort,nspin,
     &	                                    0:ndos))
      consp(1)=-1.d0
      consp(nspin)=1.d0
      deltae=(emaxdos-emindos)/ndos
      fct=evolt/2
      if(ubi==' qp') then
        call output_sym_k_orb_char(e_qp,z_qp,nbndf)
      else if(ubi=='dft'.or.ubi==' hf') then
        call output_sym_k_orb_char(e_bnd,z_bnd,nbndf)
      else
        call output_sym_k_orb_char(e_qp,z_bnd,nbndf)
      endif
c ------------------------------------------------------------------
      allocate(nd_dos(nproc_t))
      allocate(nmpi_dos(nproc_t))
      call size_shift_par(ndos+1,nproc_t,nd_dos,nmpi_dos)
      dossp=0.d0
      pdos=0.d0
      dos_k=0.d0
      if(irel/=2.and.iprn/=0) porb=0.d0
      do ispin=1,nspin
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_bnd(k,ispin)
          if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
            allocate(achain(0:n_omega,n,n))
            call a_cont_coeff(achain,ispin,ind_k,n)
          endif
          allocate(weigd(n))
          do ind_dos=1,nd_dos(me_t+1)
            iedos=nmpi_dos(me_t+1)+ind_dos-1
            om0=emindos+iedos*deltae
            e0=om0+chem_pot
c -----------------------------------------------------------------
            if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
              call sum_up_weights(weigd,om0,ispin,ind_k,n)
            else if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
              call a_cont_weights(achain,weigd,om0,ispin,ind_k,n)
            endif
c ------ Full and Partial DOS for all variants -------------------	    
            call dos(dossp(1,iedos),pdos(1,1,1,iedos),weigd,ispin,ind_k,
     &               n)
            if(irel/=2.and.iprn/=0)
     &	      call dos_orb(porb(1,1,1,iedos),weigd,ispin,ind_k,n)
c ------ K-resolved spectral function for SYM-K -------------------
            do i=1,n
              dos_k(iedos,k,ispin)=dos_k(iedos,k,ispin)+weigd(i)/wgt(k)
            enddo
          enddo   !! over ind_dos
          deallocate(weigd)
          if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') deallocate(achain)
        enddo   !! over ind_k
      enddo   !! over ispin
      deallocate(nd_dos,nmpi_dos)
      if(nproc_k/=1) then
        call DGOP(pdos,(nrel*maxb+1)*nsort*nspin*(1+ndos),'  +',comm_k)
        call DGOP(dossp,nspin*(1+ndos),'  +',comm_k)
        if(irel/=2.and.iprn/=0)
     &    call DGOP(porb,maxel_red*nsort*nspin*(1+ndos),'  +',comm_k)
        call dgop(dos_k,nspin*npnt*(1+ndos),'  +',comm_k)
      endif
      if(nproc_t/=1) then
        call dgop(dossp,nspin*(1+ndos),'  +',comm_t)
        call dgop(pdos,(nrel*maxb+1)*nsort*nspin*(1+ndos),'  +',comm_t)
        if(irel/=2.and.iprn/=0) 
     &	  call dgop(porb,maxel_red*nsort*nspin*(1+ndos),'  +',comm_t)
        call dgop(dos_k,nspin*npnt*(1+ndos),'  +',comm_t)
      endif
c --------- DOS, DOS_K, and PDOS  -------------
      if(maswrk) then
        i_len=len_trim(allfile)
c ---- Full DOS (spin UP + spin DOWN ) ---------------------------	  
        if(ubi=='dft') open(3,file=allfile(1:i_len)//'_dft.dos')
        if(ubi==' hf') open(3,file=allfile(1:i_len)//'_hf.dos')
        if(ubi=='  x') open(3,file=allfile(1:i_len)//'_x.dos')
        if(ubi==' qp') open(3,file=allfile(1:i_len)//'_qp.dos')
        if(ubi==' gw') open(3,file=allfile(1:i_len)//'_gw.dos')
        if(ubi=='psi') open(3,file=allfile(1:i_len)//'_psi.dos')
        if(ubi=='bsp') open(3,file=allfile(1:i_len)//'_bsp.dos')
        if(irel<2) then
          if(nspin==1) then
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              ds1=dossp(1,iedos)/fct/2
              write(3,'(4(f11.5,2x))')om0*fct,ds1,ds1,ds1+ds1
            enddo
          else if(nspin==2) then
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              ds1=dossp(1,iedos)/fct/2
              ds2=dossp(2,iedos)/fct/2
              write(3,'(4(f11.5,2x))')om0*fct,ds1,ds2,ds1+ds2
            enddo
          endif
        else if(irel==2) then
          do iedos=0,ndos
            om0=emindos+iedos*deltae
            ds1=dossp(1,iedos)/fct
            write(3,'(2(f11.5,2x))')om0*fct,ds1
          enddo
        endif
        close(3)
c ---- Full DOS_K (spin UP + spin DOWN ) ---------------------------
        do ispin=1,nspin
          if(ispin==1) s2='_dn'
          if(nspin==2.and.ispin==2) s2='_up'	 	  
          if(ubi=='dft') open(3,file=allfile(1:i_len)//s2//'_dft.dosk')	
          if(ubi==' hf') open(3,file=allfile(1:i_len)//s2//'_hf.dosk')
          if(ubi=='  x') open(3,file=allfile(1:i_len)//s2//'_x.dosk')
          if(ubi==' qp') open(3,file=allfile(1:i_len)//s2//'_qp.dosk')
          if(ubi==' gw') open(3,file=allfile(1:i_len)//s2//'_gw.dosk')
          if(ubi=='psi') open(3,file=allfile(1:i_len)//s2//'_psi.dosk')
          if(ubi=='bsp') open(3,file=allfile(1:i_len)//s2//'_bsp.dosk')
          write(3,'(11x,7(6x,a1,5x))')(sym0_label(i),i=1,min(n_sym,7))
          do iedos=0,ndos
            om0=emindos+iedos*deltae
            write(3,'(f11.5,7(1x,f11.5))')om0*fct,
     &       (dos_k(iedos,ksym_in_npnt(i),ispin)/fct,i=1,min(n_sym,7))
          enddo
          close(3)	 	  
          if(ubi=='dft') open(3,file=allfile(1:i_len)//s2//'_dft.doskl')
          if(ubi==' hf') open(3,file=allfile(1:i_len)//s2//'_hf.doskl')
          if(ubi=='  x') open(3,file=allfile(1:i_len)//s2//'_x.doskl')
          if(ubi==' qp') open(3,file=allfile(1:i_len)//s2//'_qp.doskl')
          if(ubi==' gw') open(3,file=allfile(1:i_len)//s2//'_gw.doskl')
          if(ubi=='psi') open(3,file=allfile(1:i_len)//s2//'_psi.doskl')
          if(ubi=='bsp') open(3,file=allfile(1:i_len)//s2//'_bsp.doskl')
          ndim=npnt/7
          if(7*ndim<npnt) ndim=ndim+1
          do ie=1,ndim
            write(3,'(11x,7(4x,i3,5x))')(i,i=(ie-1)*7+1,min(npnt,ie*7))
            do iedos=0,ndos
              om0=emindos+iedos*deltae
              write(3,'(f11.5,7(1x,f11.5))')om0*fct,
     &       (dos_k(iedos,i,ispin)/fct,i=(ie-1)*7+1,min(npnt,ie*7))
            enddo
          enddo
          close(3)
        enddo
c ---- Partial DOSs --------------------------------------------	  
        if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
          sss='___sum'
        else if(ubi==' gw'.or.ubi=='psi'.or.ubi=='bsp') then
          sss='_Chain'
        endif
        do isort=1,nsort
c -------- S P D F .... ------------------------------------	      
          if(ubi=='dft') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_dft.pdos')
          if(ubi==' hf') open(3,file=allfile(1:i_len)
     &	                 //txtel(isort)(1:2)//sss//'_hf.pdos')
          if(ubi=='  x') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_x.pdos')
          if(ubi==' gw') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_gw.pdos')
          if(ubi=='psi') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_psi.pdos')
          if(ubi=='bsp') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_bsp.pdos')
          if(ubi==' qp') open(3,file=allfile(1:i_len)
     &		                 //txtel(isort)(1:2)//sss//'_qp.pdos')
          do iedos=0,ndos
            om0=emindos+iedos*deltae
            write(3,'(9f8.3)')om0*fct,
     &			  ((pdos(li,isort,ispin,iedos)*consp(ispin)/fct,
     &                  li=1,nrel*min(3,lmb(isort))+1),ispin=1,nspin)
          enddo
          close(3)
c --- Correlated orbitals only :  S, X, Y, Z, XY ..... ------------
          if(irel/=2.and.iprn/=0) then
            do ispin=1,nspin
              if(ispin==1) s2='_DN'
              if(ispin==2) s2='_UP'
              ind0=1
              do l=0,lmb_red(isort)
                do ie=1,ntle(l,isort)
                  if(correlated(ie,l,isort)=='N') cycle
                  ndim=l+l+1
                  if(l==0) s1='_S'
                  if(l==1) s1='_P'
                  if(l==2) s1='_D'
                  if(l==3) s1='_F'
                  ss1=s1//s2
                  if(ubi=='dft') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_dft.porb')
                  if(ubi==' hf') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_hf.porb')
                  if(ubi=='  x') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_x.porb')
                  if(ubi==' gw') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_gw.porb')
                  if(ubi=='psi') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_psi.porb')
                  if(ubi=='bsp') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_bsp.porb')
                  if(ubi==' qp') open(3,file=allfile(1:i_len)
     &		             //txtel(isort)(1:2)//ss1//sss//'_qp.porb')
                  if(l==0) write(3,*)'          S '
                  if(l==1) write(3,*)'     Y   Z   X '
                  if(l==2) write(3,*)
     &               '     XY  YZ  2ZZ-XX-YY  XZ  XX-YY '
                  if(l==3) write(3,*)
     &' 3XXY  XYZ  4YZZ-XXY-YYY 2ZZZ-3XXZ-3YYZX ... '
                  do iedos=0,ndos
                    om0=emindos+iedos*deltae
                    write(3,'(9f8.3)')om0*fct,
     &			      (porb(li,isort,ispin,iedos)*nhsort(isort)/fct,
     &                  li=ind0,ind0+ndim-1)
                  enddo
                  close(3)
                enddo  !! over ie
              enddo  !! over l
            enddo  !! over ispin
          endif
        enddo   !! over isort
      endif
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      if(.not.metal) call out_band_gap(dos_k)
      deallocate(dossp,pdos,dos_k)
      if(irel/=2.and.iprn/=0) deallocate(porb)
      call timel('***** OUTPUT_DOS finished **********')
      end
      
      
      
      
      subroutine out_band_gap(dos_k)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      real*8, intent(in) :: dos_k(0:ndos,npnt,nspin)
      integer :: ispin,k,iedos,i,kv,kc
      real*8 :: emi,ema,a,b,d0,d1,deltae,om0,x
      real*8, allocatable :: lu(:),ho(:)
      if(metal) return
      emi=-10000.d0
      ema=10000.d0
      deltae=(emaxdos-emindos)/ndos
      allocate(lu(npnt))
      allocate(ho(npnt))
      lu=10.d0
      ho=-10.d0
      kv=0
      kc=0
      do ispin=1,nspin
        do k=1,npnt
          a=-10000.d0
          b=10000.d0
          if(ubi=='dft'.or.ubi==' hf') then
            do i=1,n_bnd(k,ispin)
              x=e_bnd(i,k,ispin)-chem_pot
              if(x<0) a=max(a,x)
              if(x>0) b=min(b,x)
            enddo
          else if(ubi==' qp') then
            do i=1,n_bnd(k,ispin)
              x=e_qp(i,k,ispin)-chem_pot
              if(x<0) a=max(a,x)
              if(x>0) b=min(b,x)
            enddo
          else
            d0=(dos_k(2,k,ispin)-dos_k(0,k,ispin))/2/deltae       
            do iedos=2,ndos-1
              om0=emindos+iedos*deltae
              d1=(dos_k(iedos+1,k,ispin)-dos_k(iedos-1,k,ispin))/2
     &              /deltae
              if(d0>0.d0.and.d1<=0.d0) then
                x=deltae*d0/(d0-d1)
                x=om0-deltae+x
                if(x<0) a=max(a,x)
                if(x>0) b=min(b,x)
              endif
              d0=d1
            enddo
          endif
          if(a>emi) then
            emi=a
            kv=k
          endif
          if(b<ema) then
            ema=b
            kc=k
          endif
          lu(k)=min(lu(k),b)
          ho(k)=max(ho(k),a)
        enddo
      enddo
      if(maswrk) then
        write(iun,*)' Lowest Unoccupied States (eV):'
        write(iun,'(6(i4,1x,f7.3))')(k,lu(k)*evolt/2,k=1,npnt)
        write(iun,*)' Highest Occupied States (eV):'
        write(iun,'(6(i4,1x,f7.3))')(k,ho(k)*evolt/2,k=1,npnt)
        write(iun,*)' Band Gap is ',(ema-emi)*evolt/2,' eV'
        write(iun,*)' Between k-points: ',kv,kc
      endif
      deallocate(lu,ho)
      end
