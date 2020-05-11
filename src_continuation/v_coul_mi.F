      subroutine v_coul_mi(k,vmi,npw,khbr)
c     khbr = 0 bare Coulomb interaction
c     khbr = 1 screened Coulomb interaction
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: k,npw,khbr
      complex*16, intent(out) :: vmi(ndim_tb_n_pbmt(me_t_b+1),nplw_gw)
      integer :: igg,i,ig,ie,iatom,isort,ind,i0,ind_j,j
      real*8 :: qg(3),pi8,sqpi8,pi2,pi2a,qg2,qq_2,qg3,p32,q2a
      real*8, allocatable :: mt_pw_0(:,:)
      complex*16, allocatable :: v_tmp(:,:)
c -------------------------------------------------------
      pi8=8.d0*pi
      p32=32.d0*pi**2/sqrt(amega)
      sqpi8=sqrt(pi8)
      pi2=pi+pi
      pi2a=pi2/par
      allocate(v_tmp(n_pbmt,npw))
      v_tmp=(0.d0,0.d0)
      allocate(mt_pw_0(maxlfpb,nsort))
      q2a=q2aver
      do ig=1,npw
        igg=indpw_gw(ig,k)
        do i=1,3
          qg(i)=gbs(i,igg)+pnt(i,k)
        enddo
        qg2=dot_product(qg,qg)
        if(khbr==1.or.qg2>1.d-8) then
          qg3=sqrt(qg2)
          call integral_pbmt_pw_0(mt_pw_0,qg3)
          call integral_pbmt_pw(v_tmp(1,ig),mt_pw_0,qg)
          qg2=qg2*pi2a**2
          if(khbr==0) qq_2=pi8/qg2/sqrt(amega)
          if(khbr==1) qq_2=pi8/(qg2+omega_hse**2)/sqrt(amega)
          v_tmp(:,ig)=v_tmp(:,ig)*qq_2
        endif
      enddo  !! over ig
      deallocate(mt_pw_0)
      do ind_j=1,ndim_tb_n_pbmt(me_t_b+1)
        j=nmpi_tb_n_pbmt(me_t_b+1)+ind_j
        vmi(ind_j,1:npw)=v_tmp(j,:)
      enddo
      deallocate(v_tmp)
      end
