      subroutine mbpt_spectrum
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,k,ind_k,i_len,n,i,ib,ind_tau,j,it,m,i0,i1,ii
      real*8 :: ar,br,hh,ds
      complex*16 :: w,sig,sig1
      real*8, allocatable :: eold(:,:),xold(:,:),xnew(:,:),sc(:,:),
     &                       zc(:,:),enew(:,:),sigt(:,:,:,:)
      complex*16, allocatable :: st(:,:,:),so(:,:),achain(:),t(:,:),
     &                           t1(:,:),t2(:,:),t3(:,:,:),t11(:,:)
      hh=1.d-2
      allocate(eold(nbndf_bnd,npnt))
      allocate(xold(nbndf_bnd,npnt))
      allocate(xnew(nbndf_bnd,npnt))
      allocate(sc(nbndf_bnd,npnt))
      allocate(zc(nbndf_bnd,npnt))
      allocate(enew(nbndf_bnd,npnt))
      allocate(achain(0:n_omega))
      if(maswrk) then
        i_len=len_trim(allfile)
        open(3,file=allfile(1:i_len)//'_mbpt.eig')
        write(3,'(4x,6(7x,a5))')'E_OLD','X_OLD','X_NEW','SIG_C',
     &                          '  Z  ','E_NEW'
      endif
      i0=n0_tb_ij(me_t_b+1)
      i1=i0+nd_tb_ij(me_t_b+1)
      do ispin=1,nspin
        if(maswrk) write(3,*)' Ispin = ',ispin
        eold=0.d0
        xold=0.d0
        xnew=0.d0
        sc=0.d0
        zc=0.d0
        enew=0.d0
        do ind_k=1,ndim3_k(me_k+1)
          k=n3_mpi_k(me_k+1)+ind_k
          n=n_low_bnd(k,ispin)
          allocate(st(n,2,ndim3_tau))
          allocate(so(n,0:n_omega))
          do i=1,n
            ib=ind_bands_bnd(i,k,ispin)
            eold(i,k)=e_bnd(ib,k,ispin)-chem_pot
            ii=ib*(ib+1)/2
            if(ii>i0.and.ii<=i1) then
              if(ubi_0=='dft') then
                xold(i,k)=v_xc_bb_old(ii-i0,ind_k,ispin)
              else if(ubi_0==' hf'.and.hybrid_type==0) then
                xold(i,k)=sigx_solid_k_old(ii-i0,ind_k,ispin)
              else if(ubi_0==' hf'.and.hybrid_type/=0) then
                xold(i,k)=sigx_solid_k_old(ii-i0,ind_k,ispin)
     &                   +v_xc_bb_old(ii-i0,ind_k,ispin)
              endif
              xnew(i,k)=sigx_solid_k(ii-i0,ind_k,ispin)
            endif
            do ind_tau=1,ndim3_tau
              ar=sig_c_tau(ib,ib,1,ind_tau,ind_k,ispin)
              br=sig_c_tau(ib,ib,2,ind_tau,ind_k,ispin)
              st(i,1,ind_tau)=0.5d0*(ar+br)
              st(i,2,ind_tau)=0.5d0*(br-ar)
            enddo
          enddo
          call tau_to_omega_vrt(so,st,n)
          deallocate(st)
          do i=1,n
            call pade_fermion_all(so(i,:),1,1,achain)
            w=eold(i,k)
            call pade_fermion_all_out(1,1,achain,w,sig)
            sc(i,k)=sig
            w=eold(i,k)+hh
            call pade_fermion_all_out(1,1,achain,w,sig)
            w=eold(i,k)-hh
            call pade_fermion_all_out(1,1,achain,w,sig1)
            ds=(sig-sig1)/2.d0/hh
            zc(i,k)=1.d0/(1.d0-ds)
            enew(i,k)=eold(i,k)+zc(i,k)*(xnew(i,k)+sc(i,k)-xold(i,k))
          enddo
          deallocate(so)
        enddo  !! over ind_k                !!! over k_ind
        if(nproc_k/=1) then
          call dgop(eold,nbndf_bnd*npnt,'  +',comm_k)
          call dgop(xold,nbndf_bnd*npnt,'  +',comm_k)
          call dgop(xnew,nbndf_bnd*npnt,'  +',comm_k)
          call dgop(sc,nbndf_bnd*npnt,'  +',comm_k)
          call dgop(zc,nbndf_bnd*npnt,'  +',comm_k)
          call dgop(enew,nbndf_bnd*npnt,'  +',comm_k)
        endif
        if(nproc_t*nproc_b/=1) then
          call dgop(xold,nbndf_bnd*npnt,'  +',comm_t_b)
          call dgop(xnew,nbndf_bnd*npnt,'  +',comm_t_b)
        endif
        eold=eold*evolt/2.d0
        xold=xold*evolt/2.d0
        xnew=xnew*evolt/2.d0
        sc=sc*evolt/2.d0
        enew=enew*evolt/2.d0
        if(maswrk) then
          do k=1,npnt
            n=n_low_bnd(k,ispin)
            write(3,*)' K = ',k,sym_label(k)
            do i=1,n
              ib=ind_bands_bnd(i,k,ispin)
              write(3,'(i4,6(1x,f11.5))')ib,eold(i,k),xold(i,k),
     &                                   xnew(i,k),sc(i,k),zc(i,k),
     &                                   enew(i,k)
            enddo
          enddo
        endif
      enddo  !! over ispin
      if(maswrk) close(3)
      deallocate(eold,xold,xnew,achain,sc,zc,enew)
      end
