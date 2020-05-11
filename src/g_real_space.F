      subroutine g_real_space(g_tau,g_mm,g_mi,g_im,g_ii,z_red,a_red,phf,
     &                        ispin)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: ispin
	  real*8, intent(in) :: g_tau(nbndf_bnd,nbndf_bnd,2,ndim3_tau,
     &                            npnt_c)
      complex*16, intent(in) :: z_red(nfun_red,nbndf_bnd,nqdiv_c),
     &                          a_red(nbasmpw_red,nbndf_bnd,nqdiv_c),
     &                          phf(nr_full_red,nqdiv_c)
      real*8, intent(out) :: g_mm(nrel,nfun_red,nfun_red,0:n_tau,
     &                            nqdiv_c),
     &                       g_mi(nrel,nfun_red,nrel*nr_full_red,
     &                            0:n_tau,nqdiv_c),
     &                       g_im(nrel,nrel*nr_full_red,nfun_red,
     &                            0:n_tau,nqdiv_c),
     &                       g_ii(nrel,nrel*nr_red_red,nrel*nr_full_red,
     &                            0:n_tau,nqdiv_c)
      integer :: k,ind_tau,it,nb0,i_tau,nb,k0,k1,k11,ir,ii,kf,i0,i,irl,
     &           i05,i5,nrr,nbas,nbas0,ka,ind1,ind
      complex*16 :: cc
      complex*16, allocatable :: g(:,:),tm(:,:,:),tmp(:,:),tmp1(:,:),
     &                           t1(:,:),t2(:,:)
      nrr=nrel*nr_full_red
	  nb0=max(nfun_red,nrr,nbasmpw_red)
	  allocate(g(nbndf_bnd,nbndf_bnd))
	  allocate(tm(nqdiv_c,nrs_red,nrs_red))
	  allocate(tmp(nb0,nb0))
	  allocate(tmp1(nbndf_bnd,nb0))
      allocate(t1(nr_full_red,max(nfun_red,nrr,nbndf_bnd)))
      allocate(t2(nbndf_bnd,nrr))
	  g_mm=0.d0
	  g_mi=0.d0
	  g_im=0.d0
	  g_ii=0.d0
	  do ind_tau=1,ndim3_tau
        i_tau=ndim3_tau*me_t+ind_tau-1
        do it=1,2
          if(it==2) i_tau=n_tau-i_tau
          tm=(0.d0,0.d0)
	      do k=1,npnt_c
            kf=k_a_from_c(k)
	        nb=n_low_bnd(kf,ispin)
            nbas0=nbask_red(kf)
            nbas=nbas0/nrel
c -------------------- G on K_TAU ------------------------------------
            call ferm_unpack_tau(g,g_tau(1,1,1,ind_tau,k),nb,nbndf_bnd,
     &                           nbndf_bnd,it)
            do k0=1,k_star_c(k)
	          k1=k_list_c(k0,k)
              ka=k_a_from_c(k1)
	          call zone1_number(pnt_c(1,k1),rb0_c,ndiv_c,k11)
c ------------- Common for MM and IM ---------------------------------
              call zgemm('n','c',nb,nfun_red,nb,(1.d0,0.d0),g,nbndf_bnd,
     &                   z_red(1,1,k1),nfun_red,(0.d0,0.d0),tmp1,
     &                   nbndf_bnd)
c --------------------- MM -------------------------------------------
	          call zgemm('n','n',nfun_red,nfun_red,nb,(1.d0,0.d0),
     &                   z_red(1,1,k1),nfun_red,tmp1,nbndf_bnd,
     &                   (0.d0,0.d0),tmp,nb0)
	          tm(k11,1:nfun_red,1:nfun_red)=tmp(1:nfun_red,1:nfun_red)
c --------------------- IM -------------------------------------------
	          call zgemm('n','n',nbas0,nfun_red,nb,(1.d0,0.d0),
     &                   a_red(1,1,k1),nbasmpw_red,tmp1,nbndf_bnd,
     &                   (0.d0,0.d0),tmp,nb0)
              do irl=1,nrel
                ind1=(irl-1)*nbas
                t1=(0.d0,0.d0)
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  t1(ind,1:nfun_red)=tmp(ind1+i,1:nfun_red)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),
     &                    nfun_red,t1,1)
                ind=(irl-1)*nr_full_red+nfun_red
                do i=1,nr_full_red
                  cc=phf(i,k1)/sqrt(amega)
	              tm(k11,ind+i,1:nfun_red)=cc*t1(i,1:nfun_red)
                enddo
              enddo  !! over irl
c ------------- Common for MI and II ---------------------------------
	          call zgemm('n','c',nb,nbas0,nb,(1.d0,0.d0),g,nbndf_bnd,
     &	                 a_red(1,1,k1),nbasmpw_red,(0.d0,0.d0),tmp1,
     &                   nbndf_bnd)
              do irl=1,nrel
                ind1=(irl-1)*nbas
                t1=(0.d0,0.d0)
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  t1(ind,1:nb)=tmp1(1:nb,ind1+i)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nb,t1,
     &                    -1)
                ind=(irl-1)*nr_full_red
                do i=1,nr_full_red
                  cc=conjg(phf(i,k1))/sqrt(amega)
	              t2(1:nb,ind+i)=cc*t1(i,1:nb)
                enddo
              enddo  !! over irl
c --------------------- MI -------------------------------------------
	          call zgemm('n','n',nfun_red,nrr,nb,(1.d0,0.d0),
     &                   z_red(1,1,k1),nfun_red,t2,nbndf_bnd,
     &                   (0.d0,0.d0),tmp,nb0)
	          tm(k11,1:nfun_red,nfun_red+1:nrs_red)=
     &                                            tmp(1:nfun_red,1:nrr)
c --------------------- II -------------------------------------------
	          call zgemm('n','n',nbas0,nrr,nb,(1.d0,0.d0),
     &                   a_red(1,1,k1),nbasmpw_red,t2,nbndf_bnd,
     &                   (0.d0,0.d0),tmp,nb0)
              do irl=1,nrel
                ind1=(irl-1)*nbas
                t1=(0.d0,0.d0)
                do i=1,nbas
                  i0=indgb_red(i,ka)
                  call zone1_number(gbs(:,i0),rbas,nrdiv_red,ind)
                  t1(ind,1:nrr)=tmp(ind1+i,1:nrr)
                enddo
                call fft3(nrdiv_red(1),nrdiv_red(2),nrdiv_red(3),nrr,t1,
     &                    1)
                ind=(irl-1)*nr_full_red+nfun_red
                do i=1,nr_full_red
                  cc=phf(i,k1)/sqrt(amega)
	              tm(k11,ind+i,nfun_red+1:nrs_red)=cc*t1(i,1:nrr)
                enddo
              enddo  !! over irl
	        enddo
	      enddo   !! over k
          call fft3(ndiv_c(1),ndiv_c(2),ndiv_c(3),nrs_red*nrs_red,tm,1)
          do ir=1,nqdiv_c
	        g_mm(1,:,:,i_tau,ir)=real(tm(ir,1:nfun_red,1:nfun_red))
     &                          /nqdiv_c
	        g_mi(1,:,:,i_tau,ir)=real(tm(ir,1:nfun_red,
     &                                   nfun_red+1:nrs_red))/nqdiv_c
	        g_im(1,:,:,i_tau,ir)=real(tm(ir,nfun_red+1:nrs_red,
     &                                   1:nfun_red))/nqdiv_c
            if(irel==2) then
	          g_mm(2,:,:,i_tau,ir)=imag(tm(ir,1:nfun_red,1:nfun_red))
     &                            /nqdiv_c
	          g_mi(2,:,:,i_tau,ir)=imag(tm(ir,1:nfun_red,
     &                                     nfun_red+1:nrs_red))/nqdiv_c
	          g_im(2,:,:,i_tau,ir)=imag(tm(ir,nfun_red+1:nrs_red,
     &                                     1:nfun_red))/nqdiv_c
            endif
	      enddo
          do ir=1,nqdiv_c
            do i0=1,nr_red_red
              i=r_pnt_ind_red(i0)
              do irl=1,nrel
                i05=(irl-1)*nr_red_red+i0
                i5=(irl-1)*nr_full_red+i
	            g_ii(1,i05,:,i_tau,ir)=
     &            real(tm(ir,nfun_red+i5,nfun_red+1:nrs_red))/nqdiv_c
                if(irel==2) then
	              g_ii(2,i05,:,i_tau,ir)=
     &              imag(tm(ir,nfun_red+i5,nfun_red+1:nrs_red))/nqdiv_c
                endif
              enddo
            enddo
	      enddo
        enddo  !! over it
      enddo   !! over ind_tau
	  deallocate(g,tm,tmp,tmp1,t1,t2)
	  if(nproc_t/=1) then
        ii=nqdiv_c*(n_tau+1)
	    call dgop(g_mm,nrel*nfun_red**2*ii,'  +',comm_t)
	    call dgop(g_mi,nrel**2*nfun_red*nr_full_red*ii,'  +',comm_t)
	    call dgop(g_im,nrel**2*nfun_red*nr_full_red*ii,'  +',comm_t)
        call dgop(g_ii,nrel**3*nr_full_red*nr_red_red*ii,'  +',comm_t)
      endif
c      g_ii=0.d0
      end
