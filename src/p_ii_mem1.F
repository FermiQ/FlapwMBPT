      subroutine p_ii_mem1(ind_tau,isp1,isp2,p)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau,isp1,isp2
      real*8, intent(inout) :: p(ncmpl*nd_b_pbint(me_b+1),ndim3_tn,
     &                           ndim3_k(me_k+1))
      integer :: i,ir,jrl,irl,ista,ind_ista,ind0,ista1,ista2,nnn,nsta
      integer, allocatable :: nd_b(:),n0_b(:)
      real*8, allocatable :: pg1(:,:,:)
      complex*16, allocatable :: p_c(:,:),g_c(:,:,:,:,:),gk(:,:,:),
     &                           gc(:,:,:,:,:),pgr(:,:),pg(:,:,:)
      ista1=n4_mpi_nr_red(me_k_b+1)+1
      ista2=n4_mpi_nr_red(me_k_b+1)+ndim4_nr_red(me_k_b+1)
      nnn=sum(num_r_star(ista1:ista2))
      allocate(g_c(nbasmpw,nnn,nrel,npnt,2))
      call g_ii_mem_1(1,ind_tau,g_c(1,1,1,1,1),isp1,nnn)
      call g_ii_mem_1(2,ind_tau,g_c(1,1,1,1,2),isp2,nnn)
      if(ind_tau==1) call timel('** P_II_RR_0 finished **************')
c ------------------------------------------------------------------
      allocate(pg(nplw_gw,nnn,npnt))
      ind0=0
      do ind_ista=1,ndim4_nr_red(me_k_b+1)
        ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
        nsta=num_r_star(ista)
        allocate(gc(nqdiv,nr_full,nrel,nrel,2))
        allocate(gk(nbasmpw,nrel,nqdiv))
        call g_ii_mem_2(gk,g_c(1,1,1,1,1),ista,nnn,ind0)
        call g_ii_mem_3(gk,gc(1,1,1,1,1))
        call g_ii_mem_2(gk,g_c(1,1,1,1,2),ista,nnn,ind0)
        call g_ii_mem_3(gk,gc(1,1,1,1,2))
        deallocate(gk)
c ------------------------------------------------------------------
        allocate(p_c(nqdiv,nr_full))
        p_c=(0.d0,0.d0)
        do jrl=1,nrel
          do irl=1,nrel
            do i=1,nr_full
              do ir=1,nqdiv
                p_c(ir,i)=p_c(ir,i)-gc(ir,i,irl,jrl,1)
     &                             *conjg(gc(ir,i,irl,jrl,2))
              enddo
            enddo
          enddo   !! over irl
        enddo   !! over jrl
        deallocate(gc)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,p_c,-1)
c ------------------------------------------------------------------
        allocate(pgr(nplw_gw,nqdiv))
        call p_ii_q_mem_1(p_c,pgr)
        deallocate(p_c)
        call p_ii_q_mem_2(pgr,pg,nsta,ista,nnn,ind0)
        deallocate(pgr)
        ind0=ind0+nsta
      enddo   !! over ind_ista
      deallocate(g_c)
      if(ind_tau==1) call timel('** P_II_RR_1 finished **************')
c ------------------------------------------------------------------
      allocate(nd_b(nproc_b))
      allocate(n0_b(nproc_b))
      call size_shift_par(nplw_gw,nproc_b,nd_b,n0_b)
      allocate(pg1(ncmpl*nplw_gw,nd_b(me_b+1),ndim3_k(me_k+1)))
      if(ncmpl==1) then
        call p_ii_q_mem_3r(pg,pg1,nd_b(me_b+1),n0_b(me_b+1),nnn)
      else
        call p_ii_q_mem_3(pg,pg1,nd_b(me_b+1),n0_b(me_b+1),nnn)
      endif
      deallocate(pg)
      if(ind_tau==1) call timel('** P_II_Q_MEM_3 finished ***********')
      if(ncmpl==1) then
        call p_ii_q_mem_4r(pg1,nd_b(me_b+1),n0_b(me_b+1),ind_tau,p)
      else
        call p_ii_q_mem_4(pg1,nd_b(me_b+1),n0_b(me_b+1),ind_tau,p)
      endif
      deallocate(pg1,nd_b,n0_b)
      if(ind_tau==1) call timel('** P_II_Q_MEM_4 finished ***********')
      end



      subroutine p_ii_sclp(ind_tau,isp1,isp2,pwsclp)
      use atom_mod
      use manager_mod
      use parallel_mod
      use scalapack_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: ind_tau,isp1,isp2
      real*8, intent(inout) :: pwsclp(ncmpl*nrow_pw_max,ncol_pw_max,
     &                                ndim3_tn,ndim3_k(me_k+1))
      integer :: i,ir,jrl,irl,ista,ind_ista,ind0,ista1,ista2,nnn,nsta
      integer, allocatable :: nd_b(:),n0_b(:)
      real*8, allocatable :: pg1(:,:,:)
      complex*16, allocatable :: p_c(:,:),g_c(:,:,:,:,:),gk(:,:,:),
     &                           gc(:,:,:,:,:),pgr(:,:),pg(:,:,:)
      ista1=n4_mpi_nr_red(me_k_b+1)+1
      ista2=n4_mpi_nr_red(me_k_b+1)+ndim4_nr_red(me_k_b+1)
      nnn=sum(num_r_star(ista1:ista2))
      allocate(g_c(nbasmpw,nnn,nrel,npnt,2))
      call g_ii_mem_1(1,ind_tau,g_c(1,1,1,1,1),isp1,nnn)
      call g_ii_mem_1(2,ind_tau,g_c(1,1,1,1,2),isp2,nnn)
      if(ind_tau==1) call timel('** P_II_RR_0 finished **************')
c ------------------------------------------------------------------
      allocate(pg(nplw_gw,nnn,npnt))
      ind0=0
      do ind_ista=1,ndim4_nr_red(me_k_b+1)
        ista=n4_mpi_nr_red(me_k_b+1)+ind_ista
        nsta=num_r_star(ista)
        allocate(gc(nqdiv,nr_full,nrel,nrel,2))
        allocate(gk(nbasmpw,nrel,nqdiv))
        call g_ii_mem_2(gk,g_c(1,1,1,1,1),ista,nnn,ind0)
        call g_ii_mem_3(gk,gc(1,1,1,1,1))
        call g_ii_mem_2(gk,g_c(1,1,1,1,2),ista,nnn,ind0)
        call g_ii_mem_3(gk,gc(1,1,1,1,2))
        deallocate(gk)
c ------------------------------------------------------------------
        allocate(p_c(nqdiv,nr_full))
        p_c=(0.d0,0.d0)
        do jrl=1,nrel
          do irl=1,nrel
            do i=1,nr_full
              do ir=1,nqdiv
                p_c(ir,i)=p_c(ir,i)-gc(ir,i,irl,jrl,1)
     &                             *conjg(gc(ir,i,irl,jrl,2))
              enddo
            enddo
          enddo   !! over irl
        enddo   !! over jrl
        deallocate(gc)
        call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full,p_c,-1)
c ------------------------------------------------------------------
        allocate(pgr(nplw_gw,nqdiv))
        call p_ii_q_mem_1(p_c,pgr)
        deallocate(p_c)
        call p_ii_q_mem_2(pgr,pg,nsta,ista,nnn,ind0)
        deallocate(pgr)
        ind0=ind0+nsta
      enddo   !! over ind_ista
      deallocate(g_c)
      if(ind_tau==1) call timel('** P_II_RR_1 finished **************')
c ------------------------------------------------------------------
      allocate(nd_b(nproc_b))
      allocate(n0_b(nproc_b))
      call size_shift_par(nplw_gw,nproc_b,nd_b,n0_b)
      allocate(pg1(ncmpl*nplw_gw,nd_b(me_b+1),ndim3_k(me_k+1)))
      if(ncmpl==1) then
        call p_ii_q_mem_3r(pg,pg1,nd_b(me_b+1),n0_b(me_b+1),nnn)
      else
        call p_ii_q_mem_3(pg,pg1,nd_b(me_b+1),n0_b(me_b+1),nnn)
      endif
      deallocate(pg)
      if(ind_tau==1) call timel('** P_II_Q_MEM_3 finished ***********')
      if(ncmpl==1) then
        call p_ii_q_sclp_4r(pg1,nd_b(me_b+1),n0_b(me_b+1),ind_tau,
     &                      pwsclp)
      else
        call p_ii_q_sclp_4(pg1,nd_b(me_b+1),n0_b(me_b+1),ind_tau,pwsclp)
      endif
      deallocate(pg1,nd_b,n0_b)
      if(ind_tau==1) call timel('** P_II_Q_SCLP_4 finished **********')
      end
