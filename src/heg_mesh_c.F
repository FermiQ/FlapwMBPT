      subroutine heg_mesh_c
      use atom_mod
      use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      logical :: err
      integer :: i,inf,nout1,ncosph,istar,i1,ig,j,ii,k,j0,j1,jj
      real*8 :: pi2,rout,cut,x(3),diff,pi2a,a
c     *******************************************************
      pi2=pi+pi
      gbas=0.d0
      do i=1,3
        gbas(i,i)=1.d0
      enddo
      do i=1,3
        rb0_c(:,i) = rbas(:,i)*ndiv_c(i)
        qb0_c(:,i) = gbas(:,i)/ndiv_c(i)
        r0b0_c(:,i) = rbas(:,i)/nrdiv_red(i)
        q0b0_c(:,i) = gbas(:,i)*nrdiv_red(i)
      enddo
      do i=1,3
        nr_div_c(i)=nrdiv_red(i)*ndiv_c(i)
      enddo
      nrr_div_c=nr_div_c(1)*nr_div_c(2)*nr_div_c(3)
      nqdiv_c=ndiv_c(1)*ndiv_c(2)*ndiv_c(3)
      nr_full_c=nrdiv_red(1)*nrdiv_red(2)*nrdiv_red(3)
      allocate(q_mesh_heg_c(3,nrr_div_c))
      allocate(r_mesh_heg_c(3,nrr_div_c))
      nc_line=nr_div_c(3)
      allocate(qc_line(3,nc_line))
      ii=0
      do k=0,nr_div_c(3)-1
        do j=0,nr_div_c(2)-1
          do i=0,nr_div_c(1)-1
            ii=ii+1
            q_mesh_heg_c(:,ii)=qb0_c(:,1)*i+qb0_c(:,2)*j+qb0_c(:,3)*k
            r_mesh_heg_c(:,ii)=r0b0_c(:,1)*i+r0b0_c(:,2)*j
     &                        +r0b0_c(:,3)*k
          enddo
        enddo
        if(k_line=='001') qc_line(:,k+1)=qb0_c(:,3)*k
        if(k_line=='011') qc_line(:,k+1)=(qb0_c(:,2)+qb0_c(:,3))*k
        if(k_line=='111') qc_line(:,k+1)=(qb0_c(:,1)+qb0_c(:,2)
     &                                              +qb0_c(:,3))*k
      enddo
c --------------------------------------------------------------
      allocate(gbs(3,nrr_div_c))
      allocate(lastg(0:nrr_div_c/2+1))
      allocate(rcosph(nrr_div_c/2))
      cut=qb0_c(1,1)*(nr_div_c(1)/2-1)+1.d-6
      x=0.d0
      call vgensp(cut,nrr_div_c,qb0_c,x,rout,nout1,gbs,ncosph,
     &            rcosph,lastg(1),inf,err)
      call timel('**** vgensp in heg_mesh finished ***')
      lastg(0)=0
      if(inf.eq.0.and.maswrk) then
        write(iun,1000)nout1,cutoff
        write(iun,1030)ncosph
        write(iun,*)' Number of vectors in shells (RED):'
        write(iun,'(18i4)')(lastg(i)-lastg(i-1),i=1,ncosph)
      endif
      if(inf.eq.1) then
        if(maswrk) write(iun,1010)rout,nout1
        call ending
      endif
c --------- Generating the stars ---------------------------------
      ii=3
1     allocate(indstar_c(0:ii*ncosph))
      allocate(qabs_star_c(ii*ncosph))
      call set_stars(ncosph,nstar_c,indstar_c,qabs_star_c,ii,err)
      if(err) then
        deallocate(indstar_c,qabs_star_c)
        ii=ii+1
        goto 1
      endif
      allocate(q_irr_c(3,nstar_c))
      allocate(wgt_c(nstar_c))
      do i=1,nstar_c
        q_irr_c(:,i)=gbs(:,indstar_c(i-1)+1)
        wgt_c(i)=dble(indstar_c(i)-indstar_c(i-1))/nqdiv_c
      enddo
      if(maswrk) write(iun,1050)nstar_c
c --------------------------------------------------------------
      nq_sph_c=nout1
      allocate(q_sph_c(3,nq_sph_c))
      q_sph_c=gbs(:,1:nq_sph_c)
      deallocate(gbs,lastg,rcosph)
c ----------- Correspondence to cube -----------------------------
      allocate(q_sph_in_cube_c(nq_sph_c))
      allocate(q_cube_in_sph_c(nrr_div_c))
      q_cube_in_sph_c=0
      do i=1,nq_sph_c
        call zone1_number(q_sph_c(1,i),rb0_c,nr_div_c,j)
        q_sph_in_cube_c(i)=j
        q_cube_in_sph_c(j)=i
      enddo
c --------- Symmetry relations -----------------------------------
      allocate(i_kref_c(nq_sph_c))
      allocate(k_group_c(nq_sph_c))
      allocate(ip_k_c(nstar_c,ngroup))
      do istar=1,nstar_c
        i1=indstar_c(istar-1)+1
        do i=i1,indstar_c(istar)
          i_kref_c(i)=istar
        enddo
        do ig=1,ngroup
          call rotate(q_irr_c(1,istar),q_irr_c(2,istar),
     &                q_irr_c(3,istar),x(1),x(2),x(3),u(2,ig),1)
          do i=i1,indstar_c(istar)
            diff=maxval(abs(q_sph_c(:,i)-x))
            if(diff<1.d-8) then
              ip_k_c(istar,ig)=i
              k_group_c(i)=ig
            endif
          enddo
        enddo
      enddo
c ----------- Correspondence to line -----------------------------
      allocate(kline_in_npnt_c(nc_line))
      allocate(k_npnt_in_line_c(nstar_c))
      k_npnt_in_line_c=0
      kline_in_npnt_c=0
      jj=nc_line
      do i=1,nc_line
        call zone1_number(qc_line(1,i),rb0_c,nr_div_c,j)
        j1=q_cube_in_sph_c(j)
        if(j1==0) then
          jj=i-1
          goto 2
        endif
        j0=i_kref_c(j1)
        kline_in_npnt_c(i)=j0
        k_npnt_in_line_c(j0)=i
      enddo
2     nc_line=jj
      allocate(qabs_line_c(nc_line))
      do i=1,nc_line
        j0=kline_in_npnt_c(i)
        qabs_line_c(i)=qabs_star_c(j0)
      enddo
      allocate(ndim_k_red(nproc_k))
      allocate(n_mpi_k_red(nproc_k))
      call size_shift_par(nstar_c,nproc_k,ndim_k_red,n_mpi_k_red)
c ---------- Closest to K_fermi K-point ---------------------------
      pi2a=pi2/par
      n_k_fermi_c=1
      diff=k_fermi/pi2a
      do i=2,nstar_c
        a=abs(k_fermi/pi2a-qabs_star_c(i))
        if(a<diff) then
          diff=a
          n_k_fermi_c=i
        endif
      enddo
      n_kl_fermi_c=1
      diff=k_fermi/pi2a
      do i=2,nc_line
        a=abs(k_fermi/pi2a-qabs_line_c(i))
        if(a<diff) then
          diff=a
          n_kl_fermi_c=i
        endif
      enddo
      if(maswrk) write(iun,*)' N_K_FERMI_C, N_KL_FERMI_C:',n_k_fermi_c,
     &                       n_kl_fermi_c
      call timel('*** HEG_MESH_C finished ************')
c -----------------------------------------------------------------
1000  format(1x,'All ',i6,' vectors in sphere ',f9.4,' generated OK')
1010  format(1x,'WARNING from GETPLW : In sphere ',f9.4,' generated ',
     &       i6,' vectors')
1030  format(1x,'Total # of co-spheres for RO and V found is ',i7) 
1050  format(1x,'Total # of stars of q-vectors found is ',i7)     
      end
