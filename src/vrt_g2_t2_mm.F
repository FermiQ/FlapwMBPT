      subroutine vrt_g2_t2_mm(ispin,iq,krs,z_red,trr,lambda_tau,ch)
      use atom_mod
      use etot_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      character*2, intent(in) :: ch
      integer, intent(in) :: iq,ispin
      real*8, intent(in) :: trr(n_lem2_red,n_lem2_red,ndimc_kk(me_k+1))
      complex*16, intent(in) :: krs(nfun_red,nfun_red,ndimc_kk(me_k+1)),
     &                          z_red(nfun_red,nbndf_bnd,nqdiv_c,nspin)
      complex*16, intent(inout) :: lambda_tau(nbndf_bnd,nbndf_bnd,
     &                                        ndimc_kk(me_k+1))
      integer :: k,i,ij,jatom,iatom,jsort,isort,ind,ind2,ind_r,ir,jnd,
     &           jnd2,j,kl,n,n1,l,nn,nn1
      complex*16, allocatable :: vx(:,:,:)
      allocate(vx(nfun_red,nfun_red,nqdiv_c))
      vx=(0.d0,0.d0)
      do jatom=1,natom
        jsort=is(jatom)
        n1=lfunm_red(jsort)
        nn1=n1*n1
        jnd=io_lem_red(jatom)-1
        jnd2=io_lem2_red(jatom)-1
        do iatom=1,natom
          isort=is(iatom)
          n=lfunm_red(isort)
          nn=n*n
          ind=io_lem_red(iatom)-1
          ind2=io_lem2_red(iatom)-1
          do ind_r=1,ndimc_kk(me_k+1)
            ir=nc_mpi_kk(me_k+1)+ind_r
            kl=0
            do l=1,n1
              do k=1,n1
                kl=kl+1
                ij=0
                do j=1,n
                  do i=1,n
                    ij=ij+1
                    vx(ind+j,jnd+l,ir)=vx(ind+j,jnd+l,ir)
     &                                +trr(ind2+ij,jnd2+kl,ind_r)
     &                                *krs(ind+i,jnd+k,ind_r)
                  enddo
                enddo
              enddo
            enddo
          enddo  !! over ind_r
        enddo  !! over iatom
      enddo  !! over jatom
      if(nproc_k/=1) call dgop(vx,2*nfun_red**2*nqdiv_c,'  +',comm_k)
      call vertex_rr_to_k_all(vx,nfun_red**2)
      call vertex_sf_to_bnd_mm(vx,lambda_tau,iq,z_red,ispin,ch)
      deallocate(vx)
      end
